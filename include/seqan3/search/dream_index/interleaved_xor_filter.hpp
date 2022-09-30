/*!\file
 * \author Jens-Uwe Ulrich <jens-uwe.ulrich AT hpi.de>
 * \brief Provides seqan3::interleaved_xor_filter.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/bit>

#include <sdsl/bit_vectors.hpp>
#include <bitset>

#include <seqan3/core/concept/cereal.hpp>

namespace seqan3
{


struct t2val {
  uint64_t t2 = 0;
  uint64_t t2count = 0;
};

typedef struct t2val t2val_t;
const int blockShift = 18;

/*!\brief The IXF binning directory. A data structure that efficiently answers set-membership queries for multiple bins.
 * \ingroup search_dream_index
 * \implements seqan3::cerealisable
 *
 * \details
 *
 * ### Binning Directory
 *
 * A binning directory is a data structure that can be used to determine set membership for elements.
 * For example, a common use case is dividing a database into a fixed number (e.g. 1024) bins by some means
 * of clustering (e.g. taxonomic binning or k-mer similarity clustering for genomic sequences).
 * For a query, the binning directory can now answer in which bins the query (probably) occurs.
 * In SeqAn we provide the Interleaved Bloom Filter (IBF) that can answer these queries efficiently.
 *
 * ### Interleaved XOR Filter (IXF)
 *
 * TODO: add description
 * The implementation of the XOR filter construction and query is based on the fastfilter_cpp github repository
 * https://github.com/FastFilter/fastfilter_cpp
 *
 * The Interleaved XOR Filter now applies the concept of a XOR Filter to multiple sets and provides a *global*
 * data structure to determine set membership of a query in `b` data sets/bins.
 * Conceptually, a XOR Filter is created for each bin using the same fixed length and fixed hash functions for each
 * filter. The resulting `b` XOR Filters are then interleaved such that the `i`'th l-bits if each XOR Filter are
 * adjacent to each other, with l being 8 or 16 . For l=2: 
 * ```
 *  XOR Filter 0         XOR Filter 1        XOR Filter 2        XOR Filter 3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * Where `x.y` denotes the `y`'th bit of the `x`'th XOR Filter.
 * ```
 * Interleaved XOR Filter
 * |0.0|0.1|1.0|1.1|2.0|2.1|3.0|3.1|0.2|0.3|1.2|1.3|2.2|2.3|3.2|3.3|
 * ```
 * A query can now be searched in all `b` bins by computing the `h` hash functions, retrieving the `h` sub-bitvectors of
 * length `b * l` starting at the positions indicated by the hash functions. The bitwise XOR of these sub-bitvectors yields
 * the binningvector, a bitvector of length `b * l` where the `i*l`'th bits indicate set membership in the `i`'th bin iff all 
 * `l` bits of `i` are 0.
 *
 * ### Querying
 * To query the Interleaved XOR Filter for a value, call seqan3::interleaved_xor_filter::membership_agent() and use
 * the returned seqan3::interleaved_xor_filter::membership_agent.
 *
 * To count the occurrences of a range of values in the Interleaved XOR Filter, call
 * seqan3::interleaved_xor_filter::counting_agent() and use
 * the returned seqan3::interleaved_xor_filter::counting_agent_type.
 *
 *
 * ### Thread safety
 *
 * The Interleaved XOR Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time). Furthermore, XOR Filters are immutable data structures
 * by itself.
 *
 * 
 */
template <typename FingerprintType = uint8_t>
//!\cond
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
//!\endcond
class interleaved_xor_filter
{
private:
    
    using data_type = sdsl::int_vector<>;

    //!\brief The number of bins specified by the user.
    size_t bins{};
    //!\brief The number of bins stored in the IXF (next multiple of 64 of `bins`).
    //size_t technical_bins{};
    //!\brief The size of each bin in bits.
    size_t bin_size_{};
    //!\brief number of elements that can be stored in each bin.
    size_t max_bin_elements{};
    //!\brief The number of 64-bit integers needed to store `bins` many bits (e.g. `bins = 50` -> `bin_words = 1`).
    size_t bin_words{};
    //!\brief The length of each xor_filter block (filter_size/3)
    size_t block_length{};
    //!\brief The int vector of fingerprints
    data_type data{};
    //!\brief number of bits used to store one hashed item in the XOR filter
    size_t ftype{};
    //!\brief seed for  hashing
    size_t seed{13572355802537770549ULL};
    //!\brief number of bins to query in parallel
    size_t bins_per_batch{};

    /*!\brief Utilizes hashing of a 64-bit key.

     * \param   64-bit key to hash
     * \returns a 64-bit hash value for the given key
     * 
     */
    inline constexpr uint64_t murmur64(uint64_t h) const
    {
        h += seed;
        h ^= h >> 33;
        h *= UINT64_C(0xff51afd7ed558ccd);
        h ^= h >> 33;
        h *= UINT64_C(0xc4ceb9fe1a85ec53);
        h ^= h >> 33;
        return h;
    }

    /*!\brief Rotates the the bits of 64-bit hash value to the left
     * \param   64-bit hash value
     * \param   number of bits to rotate to the left
     * 
     * \returns the rotated 64-bit hash value
     * 
     */
    inline constexpr uint64_t rotl64(uint64_t n, unsigned int c) const
    {
        // assumes width is a power of 2
        const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
        // assert ( (c<=mask) &&"rotate by type width or more");
        c &= mask;
        return (n << c) | ( n >> ((-c) & mask));
    }

    /*!\brief   Fingerprint funtion of the underlying IXF
     * \param   64-bit hash value
     *  
     * \returns an unsigned 8-bit or 16-bit unsigned integer fingerprint for 
     * storing hash values in the interleaved XOR filter
     * 
     */
    inline constexpr FingerprintType fingerprint(const uint64_t hash) const 
    {
        FingerprintType h = (FingerprintType) hash ^ (hash >> 32);
        if (h == 0)
            h = (FingerprintType) hash ^ (hash >> 16);
        if (h == 0)
            h = (FingerprintType) hash ^ (hash >> 8);
        return h;
    }

    /*!\brief   Multiply hash value with blockLength and divide by 2^32
     * \param hash 32-bit hash value
     * \param n block length used for IXF construction
     *  
     * corresponds to (hash * n) / 2^32 
     *  but is 4 times faster than ordinary division
     *  more efficient use of modulo operation
     * 
     */
    inline constexpr uint32_t reduce(uint32_t hash, uint32_t n) const
    {
        // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
        return (uint32_t) (((uint64_t) hash * n) >> 32);
    }

    /*!\brief Calculate a hash from a given 64-bit hash value
     * \param hash 64-bit hash value
     * \param index number of the used hash function (1, 2 or 3)
     * \param blockLength used for IXF construction
     *  
     */
    inline constexpr size_t getHashFromHash(uint64_t hash, int index, int blockLength) const
    {
        uint32_t r = rotl64(hash, index * 21);
        return (size_t) reduce(r, blockLength) + ((size_t)index) * ((size_t)blockLength);
    }

    
    /*!\brief Count number of occurences of each index in a bin
     * \param tmp pointer to an sdsl::int_vector<> data structure
     * \param b Block number in which we count
     * \param len corresponds to the number of keys that match to the index
     * \param t2vals data structure to store counts for each index
     *  
     */
    void applyBlock(uint64_t* tmp, int b, int len, t2val_t* t2vals)
    {
        for (int i = 0; i < len; i += 2) {
            uint64_t x = tmp[(b << blockShift) + i];
            int index = (int) tmp[(b << blockShift) + i + 1];
            t2vals[index].t2count++;
            t2vals[index].t2 ^= x;
        }
    }

    /*!\brief Remove index from index occurence array and add index to single entry set
     * \param tmp pointer to an sdsl::int_vector<> data structure
     * \param b Block number in which we count
     * \param len corresponds to the number of keys that match to the index
     * \param t2vals data structure to store counts for each index
     * \param alone reference to an array of indexes that occure only once in t2vals
     * \param alonePos number of indexes that occur only once in t2vals
     *  
     * remove all indexes in tmp array from index occurence array and remove corresponding hashes by XOR-ing with the hash value
     *  add index to single entry set if an index now occurs only once in the remaining set of indexes
     */
    int applyBlock2(uint64_t* tmp, int b, int len, t2val_t*  t2vals, sdsl::int_vector<>& alone, int alonePos) 
    {
        for (int i = 0; i < len; i += 2) {
            uint64_t hash = tmp[(b << blockShift) + i];
            int index = (int) tmp[(b << blockShift) + i + 1];
            int oldCount = t2vals[index].t2count;
            if (oldCount >= 1) {
                int newCount = oldCount - 1;
                t2vals[index].t2count = newCount;
                if (newCount == 1) {
                    alone[alonePos++] = index;
                }
                t2vals[index].t2 ^= hash;
            }
        }
        return alonePos;
    }

    /*!\brief Find all index positions to which only one hashed key maps
     * \param elements list of keys to store in the XOR filter
     * \param t2vals data structure to store counts for each index
     * \param alone_positions reference to an array of indexes that occure only once
     *  
     */
    int find_alone_positions(std::vector<size_t>& elements, std::vector<t2val_t>& t2vals, sdsl::int_vector<>& alone_positions)
    {
        std::fill(t2vals.begin(), t2vals.end(), t2val_t{ 0,0 });
        // number of elements / 2^18 => if more than 2^18 elements, we need 2 blocks
        int blocks = 1 + ((3 * block_length) >> blockShift);
        sdsl::int_vector<> tmp = sdsl::int_vector(blocks << blockShift, 0, 64);       
        sdsl::int_vector<> tmpc = sdsl::int_vector(blocks,0);
        
        for(size_t k : elements) {
            uint64_t hash = murmur64(k);
            for (int hi = 0; hi < 3; hi++) {
                int index = getHashFromHash(hash, hi, block_length);
                int b = index >> blockShift;
                int i2 = tmpc[b];
                tmp[(b << blockShift) + i2] = hash;
                tmp[(b << blockShift) + i2 + 1] = index;
                tmpc[b] += 2;
                if (i2 + 2 == (1 << blockShift)) {
                    applyBlock(tmp.data(), b, i2 + 2, t2vals.data());
                    tmpc[b] = 0;
                }
            }
        }
       
        // count occurences of index positions for all computed hash values
        for (int b = 0; b < blocks; b++) {
            applyBlock(tmp.data(), b, tmpc[b], t2vals.data());
        }

        // pick only index positions where only one unique hash value points to => those are our start positions
        int alonePos = 0;
        for (size_t i = 0; i < bin_size_; i++) 
        {
            if (t2vals[i].t2count == 1) {
                alone_positions[alonePos++] = i;
            }
        }
        return alonePos;
    }

    /*!\brief Determine the order in which keys are added to the filter
     * \param reverse_order resulting stack of keys to insert in the filter
     * \param reverse_h stack of hash function numbers
     * \param t2vals data structure to store counts for each index
     * \param alone_positions reference to an array of indexes that occure only once
     * \param size number of keys to insert
     * \param alone_pos number of indexes that occur only once in t2vals
     * 
     * \returns the number of keys in the stack
     *  
     */
    size_t fill_stack(sdsl::int_vector<>& reverse_order, 
                      sdsl::int_vector<>& reverse_h, 
                      std::vector<t2val_t>& t2vals, 
                      sdsl::int_vector<>& alone_positions,
                      size_t size, 
                      int alone_pos)
    {
        int blocks = 1 + ((3 * block_length) >> blockShift);
        sdsl::int_vector<> tmp = sdsl::int_vector(blocks << blockShift, 0, 64);
        sdsl::int_vector<> tmpc = sdsl::int_vector(blocks, 0);
        size_t reverse_order_pos = 0;
        int best_block = 0;
       
        while (reverse_order_pos < size)
        {
            if (alone_pos == 0) 
            {
                // we need to apply blocks until we have an entry that is alone
                // (that is, until alonePos > 0)
                // so, find a large block (the larger the better)
                // but don't need to search very long
                // start searching where we stopped the last time
                // (to make it more even)
                
                for (int i = 0, b = best_block, best = -1; i < blocks; i++)
                {
                    if (b >= blocks)
                    {
                        b = 0;
                    }
                    if (tmpc[b] > best) 
                    {
                        best = tmpc[b];
                        best_block = b;
                        if (best > (1 << (blockShift - 1))) 
                        {
                            // sufficiently large: stop
                            break;
                        }
                    }
                }
                
                if (tmpc[best_block] > 0) {
                    alone_pos = applyBlock2(tmp.data(), best_block, tmpc[best_block], t2vals.data(), alone_positions, alone_pos);
                    tmpc[best_block] = 0;
                }
                // applying a block may not actually result in a new entry that is alone
                if (alone_pos == 0) {
                    for (int b = 0; b < blocks && alone_pos == 0; b++) {
                        if (tmpc[b] > 0) {
                            alone_pos = applyBlock2(tmp.data(), b, tmpc[b], t2vals.data(), alone_positions, alone_pos);
                            tmpc[b] = 0;
                        }
                    }
                }
            }
            if (alone_pos == 0) {
                break;
            }
            
            int i = alone_positions[--alone_pos];
            int b = i >> blockShift;
            if (tmpc[b] > 0) {
                alone_pos = applyBlock2(tmp.data(), b, tmpc[b], t2vals.data(), alone_positions, alone_pos);
                tmpc[b] = 0;
            }
            
            uint8_t found = -1;
            if (t2vals[i].t2count == 0) {
                continue;
            }
            long hash = t2vals[i].t2;
            
            for (int hi = 0; hi < 3; hi++) {
                int h = getHashFromHash(hash, hi, block_length);
                if (h == i) {
                    found = (uint8_t) hi;
                    t2vals[i].t2count = 0;
                } else {
                    int b = h >> blockShift;
                    int i2 = tmpc[b];
                    tmp[(b << blockShift) + i2] = hash;
                    tmp[(b << blockShift) + i2 + 1] = h;
                    tmpc[b] +=  2;
                    if (tmpc[b] >= 1 << blockShift) {
                        alone_pos = applyBlock2(tmp.data(), b, tmpc[b], t2vals.data(), alone_positions, alone_pos);
                        tmpc[b] = 0;
                    }
                }
            }
            
            reverse_order[reverse_order_pos] = hash;
            reverse_h[reverse_order_pos] = found;
            reverse_order_pos++;
        }
        return reverse_order_pos;
    }

    /*!\brief Adds a 8 or 16 bit value to the xor filter for each hashed key
     * \param reverse_order stack of keys to insert in the filter
     * \param reverse_h stack of hash function numbers
     * \param bin number of the interleaved XOR filter
     * 
     *  
     */
    void fill_filter(sdsl::int_vector<>& reverse_order, sdsl::int_vector<>& reverse_h, uint bin)
    {
        
        for (int i = reverse_order.size() - 1; i >= 0; i--) 
        {
            // the hash of the key we insert next
            uint64_t hash = reverse_order[i];
            int found = reverse_h[i];
            // which entry in the table we can change
            int change = -1;
            // we set table[change] to the fingerprint of the key,
            // unless the other two entries are already occupied
            FingerprintType xor2 = fingerprint(hash);
            for (int hi = 0; hi < 3; hi++) 
            {
                size_t h = getHashFromHash(hash, hi, block_length);
                if (found == hi) 
                {
                    change = h;
                } 
                else 
                {
                    uint64_t idx = bins * h;
                    idx += bin;
                    xor2 ^= data[idx];

                }
            }
            uint64_t idx = bins * change;
            idx += bin;
           
            data[idx] = xor2;
        }
    }
    
    /*!\brief Adds all keys of all bins in one go
     * \param elements two-dimensional list of bins and their keys
     *  
     */
    void add_elements(std::vector<std::vector<size_t>>& elements)
    {
        // stack sigma
        // order in which elements will be inserted into their corresponding bins
        std::vector<sdsl::int_vector<>> reverse_orders;
        // order in which hash seeds are used for element insertion
        std::vector<sdsl::int_vector<>> reverse_hs;
        size_t reverse_order_pos;
        // repeat until all xor filters can be build with same seed
        while (true)
        {
            // sets the same seed for all xor filters
            int i = 0;
            bool success = true;
            reverse_hs.clear();
            reverse_order_pos = 0;
            reverse_orders.clear();
            for (std::vector<size_t> vec : elements)
            {
                std::vector<t2val_t> t2vals_vec(bin_size_);
                sdsl::int_vector<> alone_positions = sdsl::int_vector(bin_size_);
                int alone_position_nr = find_alone_positions(vec, t2vals_vec, alone_positions);
               
                sdsl::int_vector<> rev_order_i = sdsl::int_vector(vec.size(),0,64);
                sdsl::int_vector<> reverse_hi = sdsl::int_vector(vec.size(),0,8);
                reverse_order_pos = fill_stack(std::ref(rev_order_i), std::ref(reverse_hi), t2vals_vec, alone_positions, vec.size(), alone_position_nr);
                reverse_orders.emplace_back(std::move(rev_order_i));
                reverse_hs.emplace_back(std::move(reverse_hi));
                if (reverse_order_pos != vec.size())
                {
                    success = false;
                    set_seed();
                    break;
                }

                i++;
            }
            if (success)
                break;
        }
        
        for (size_t i = 0; i < elements.size(); ++i)
        {
            fill_filter(reverse_orders[i], reverse_hs[i], i);
        }
        
    }

public:

    class membership_agent; // documented upon definition below
    template <std::integral value_t>
    class counting_agent_type; // documented upon definition below

    /*!\name Constructors, destructor and assignment
     */
    interleaved_xor_filter() = default; //!< Defaulted.
    interleaved_xor_filter(interleaved_xor_filter const &) = default; //!< Defaulted.
    interleaved_xor_filter & operator=(interleaved_xor_filter const &) = default; //!< Defaulted.
    interleaved_xor_filter(interleaved_xor_filter &&) = default; //!< Defaulted.
    interleaved_xor_filter & operator=(interleaved_xor_filter &&) = default; //!< Defaulted.
    ~interleaved_xor_filter() = default; //!< Defaulted.

    /*!\brief Construct an Interleaved XOR Filter from given sets of keys
     * \param elements two-dimensional list of bins and keys to store in the IXF
     *
     * \attention This constructor can only be used to construct Interleaved XOR Filters from 
     * the given sets of keys. Further adding of bins/keys will corrupt the data structure.
     *
     * \details
     *
     * ### Example
     * 
     */
    interleaved_xor_filter(std::vector<std::vector<size_t>>& elements)
    {
        bins = elements.size();
        for (std::vector<size_t> v : elements)
        {
            if (v.size() > max_bin_elements)
                max_bin_elements = v.size();
        }
        bin_size_ = 32 + 1.23 * max_bin_elements;
        block_length = bin_size_ / 3;
        ftype = CHAR_BIT * sizeof(FingerprintType);
        bins_per_batch = 64/ftype;

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (bin_size_ == 0)
            throw std::logic_error{"The size of a bin must be > 0."};


        if (ftype == 8)
        {
            bin_words = (bins + 7) >> 3; // = ceil(bins/8)
        }
        else
        {
            bin_words = (bins + 15) >> 2; // = ceil(bins/4)
        }

        
        data = sdsl::int_vector<>(bins * bin_size_, 0, ftype);
    
        add_elements(elements);

        
    }

    /*!\brief Construct an Interleaved XOR Filter object for later adding of bins
     * \param bins_ The number of bins.
     * \param size  maximum number of elements to store in a single filter
     *
     * \attention This constructor should be used if large numbers of bins and keys shall be
     * stored in an interleaved XOR filter. It creates the necessary data structure based
     * on the given number of bins and elements, which can be added later on by utilizing
     * add_bin_elements() function
     *
     * \details
     *
     * ### Example
     * 
     */
    interleaved_xor_filter(size_t bins_, size_t max_bin_elements)
    {
        this->bins = bins_;
        bin_size_ = 32 + 1.23 * max_bin_elements;
        block_length = bin_size_ / 3;
        ftype = CHAR_BIT * sizeof(FingerprintType);
        bins_per_batch = 64/ftype;

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (max_bin_elements == 0)
            throw std::logic_error{"The number of elements to store must be > 0."};


        if (ftype == 8)
        {
            bin_words = (bins + 7) >> 3; // = ceil(bins/8)
        }
        else
        {
            bin_words = (bins + 15) >> 2; // = ceil(bins/4)
        }

        data = sdsl::int_vector<>(bins * bin_size_, 0, ftype);
    }

    /*!\brief add all keys of a vector/list to the given bin number
     * \param bin       bins number of the interleaved XOR filter
     * \param elements  elements to store in the given bin
     * \returns `true` if elements could be added successfully, `false` otherwise.
     * 
     * \attention If this method returns `false`, you have to clear the filter (@see clear()),
     * reset the seed (@see set_seed()) and repeat adding all elements of all bins. 
     *
     * \details
     *
     * ### Example
     * 
     *  seqan3::interleaved_xor_filter<> ixf(100, 2000);
	 *  std::vector<uint64_t> elems{};
	 *  while (true)
	 *  {
     *      bool success = true;
	 *	    for (int e = 0; e < 100 ; ++e)
	 *	    {
	 *		    std::vector<uint64_t> tmp{};
	 *		    for (uint64_t i = 0; i < 2000; ++i)
	 *		    {
	 * 			    uint64_t key = (e*2000) + i;
	 *			    tmp.emplace_back(key);
	 *		    }
	 *		    success = ixf.add_bin_elements(e, tmp);
	 *		    if (!success)
	 *		    {
	 *			    ixf.clear();
	 *			    ixf.set_seed();
	 *			    break;
	 *		    }
	 *		    if (e == 2)
	 *			    elems=std::move(tmp);
	 *	    }
     *
	 *	    if (success)
	 *		    break;
	 *  }
     * 
     */
    bool add_bin_elements(size_t bin, std::vector<size_t>& elements)
    {
       
        size_t reverse_order_pos = 0;
        std::vector<t2val_t> t2vals_vec(bin_size_);
        sdsl::int_vector<> alone_positions = sdsl::int_vector(bin_size_);
        int alone_position_nr = find_alone_positions(elements, t2vals_vec, alone_positions);
        sdsl::int_vector<> rev_order = sdsl::int_vector(elements.size(),0,64);
        sdsl::int_vector<> reverse_h = sdsl::int_vector(elements.size(),0,8);
        reverse_order_pos = fill_stack(std::ref(rev_order), std::ref(reverse_h), t2vals_vec, alone_positions, elements.size(), alone_position_nr);
        if (reverse_order_pos != elements.size())
        {
            return false;
        }

        fill_filter(rev_order, reverse_h, bin);

        return true;
    }
   
    /*!\brief Clears the entire IXF.
     * \
     * \attention This function only sets the entries in the filter to 0. It will not
     * resize the filter or reset the number of bins.
     *
     * \details
     *
     * ### Example
     *
     */
    void clear()
    {
        for (size_t idx = 0; idx < data.size(); ++idx)
            data[idx] = 0;
    }

    /*!\brief Sets a new seed for the hashing function
     * \
     * \attention Setting a new seed requires clearing the entire filter and repeat
     * adding all elements of all bins.
     *
     * \details
     *
     * ### Example
     *
     */
    void set_seed()
    {
        ::std::random_device random;
        seed = random();
        seed <<= 32;
        seed |= random();
    }

    /*!\name Lookup
     * \{
     */
    /*!\brief Returns a seqan3::interleaved_xor_filter::membership_agent to be used for lookup.
     * `seqan3::interleaved_xor_filter::membership_agent`s constructed for this Interleaved XOR Filter.
     *
     * \details
     *
     * ### Example
     *
     */
    membership_agent membership_agent() const
    {
        return typename interleaved_xor_filter<FingerprintType>::membership_agent{*this};
    }

    /*!\brief Returns a seqan3::interleaved_xor_filter::counting_agent_type to be used for counting.
     * `seqan3::interleaved_xor_filter::counting_agent_type`s constructed for this Interleaved XOR Filter.
     *
     * \details
     *
     * ### Example
     *
     */
    template <typename value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
    

    /*!\brief Returns the number of bins that the Interleaved XOR Filter manages.
     * \returns The number of bins.
     */
    size_t bin_count() const noexcept
    {
        return bins;
    }

    /*!\brief Returns the size of a single bin that the Interleaved XOR Filter manages.
     * \returns The size in bits of a single bin.
     */
    size_t bin_size() const noexcept
    {
        return bin_size_;
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    size_t bit_size() const noexcept
    {
        return data.bit_size();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Test for equality.
     * \param[in] lhs A `seqan3::interleaved_xor_filter`.
     * \param[in] rhs `seqan3::interleaved_xor_filter` to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    friend bool operator==(interleaved_xor_filter const & lhs, interleaved_xor_filter const & rhs) noexcept
    {
        return std::tie(lhs.bins, lhs.bins_per_batch, lhs.bin_size_, lhs.ftype, lhs.bin_words, lhs.block_length,
                        lhs.seed, lhs.data) ==
               std::tie(rhs.bins, rhs.bins_per_batch, rhs.bin_size_, rhs.ftype, rhs.bin_words, rhs.block_length,
                        rhs.seed, rhs.data);
    }

    /*!\brief Test for inequality.
     * \param[in] lhs A `seqan3::interleaved_xor_filter`.
     * \param[in] rhs `seqan3::interleaved_xor_filter` to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    friend bool operator!=(interleaved_xor_filter const & lhs, interleaved_xor_filter const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(bins);
        archive(bin_size_);
        archive(bins_per_batch);
        archive(bin_words);
        archive(ftype);
        if (ftype != (CHAR_BIT * sizeof(FingerprintType)))
        {
            throw std::logic_error{"The interleaved XOR filter was built with a fingerprint size of " + std::to_string(ftype) +
                                   " but it is being read into an interleaved XOR filter with fingerprint of size " +
                                   std::to_string(CHAR_BIT * sizeof(FingerprintType)) + "."};
        }
        archive(block_length);
        archive(seed);
        archive(data);
    }
    //!\endcond
};

/*!\brief Manages membership queries for the seqan3::interleaved_xor_filter.
 * 
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/membership_agent_construction.cpp
 */
template <typename FingerprintType>
class interleaved_xor_filter<FingerprintType>::membership_agent
{
private:
    //!\brief The type of the augmented seqan3::interleaved_xor_filter.
    using ixf_t = interleaved_xor_filter<FingerprintType>;

    //!\brief A pointer to the augmented seqan3::interleaved_xor_filter.
    ixf_t const * ixf_ptr{nullptr};

public:
    class binning_bitvector;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default; //!< Defaulted.
    membership_agent(membership_agent const &) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default; //!< Defaulted.
    ~membership_agent() = default; //!< Defaulted.

    /*!\brief Construct a membership_agent from a seqan3::interleaved_xor_filter.
     * \private
     * \param ixf The seqan3::interleaved_xor_filter.
     */
    explicit membership_agent(ixf_t const & ixf) :
        ixf_ptr(std::addressof(ixf)), result_buffer(ixf.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    binning_bitvector result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of a given value.
     * \param[in] value The raw value to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/membership_agent_bulk_contains.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_xor_filter::membership_agent for each thread.
     */

     [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) & noexcept
    {
        assert(ixf_ptr != nullptr);
        assert(result_buffer.size() == ixf_ptr->bin_count());

        size_t bins = ixf_ptr->bin_count();
        uint8_t ftype = ixf_ptr->ftype;
        uint8_t bins_per_batch = ixf_ptr->bins_per_batch;
        uint64_t hash = ixf_ptr->murmur64(value);
        FingerprintType f = ixf_ptr->fingerprint(hash);
        uint32_t r0 = (uint32_t) hash;
        uint32_t r1 = (uint32_t) ixf_ptr->rotl64(hash, 21);
        uint32_t r2 = (uint32_t) ixf_ptr->rotl64(hash, 42);
        uint32_t h0 = ixf_ptr->reduce(r0, ixf_ptr->block_length);
        uint32_t h1 = ixf_ptr->reduce(r1, ixf_ptr->block_length) + ixf_ptr->block_length;
        uint32_t h2 = ixf_ptr->reduce(r2, ixf_ptr->block_length) + 2 * ixf_ptr->block_length;

        h0 = (h0*bins) * ftype;
        h1 = (h1*bins) * ftype;
        h2 = (h2*bins) * ftype;

        // concatenate (64/ftype) the fingerprint
        uint64_t fc64 = f;
        for (uint8_t b = 1; b < bins_per_batch ; ++b)
        {
            fc64 = (fc64 << ftype) | f;
        }

        for (size_t batch = 0; batch < ixf_ptr->bin_words; ++batch)
        {
            size_t batch_start = batch * 64;
            uint64_t v = fc64 ^ ixf_ptr->data.get_int(h0 + batch_start, 64) 
                              ^ ixf_ptr->data.get_int(h1 + batch_start, 64) 
                              ^ ixf_ptr->data.get_int(h2 + batch_start, 64);

            size_t used_bins = batch * bins_per_batch;
            uint8_t bits = 0;
            for (size_t bin = 0; bin < bins_per_batch; ++bin)
            {
                if (used_bins + bin == result_buffer.size())
                    break;

                uint64_t tmp = v << ((bins_per_batch - (bin+1)) * ftype ); 
                uint8_t tmpb = std::bitset<8>(tmp >> (64-ftype)).none() << bin;
                bits |= tmpb;
                   
            }
            result_buffer.data.set_int(used_bins, bits, ftype);
        }

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    //!\}
     [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) && noexcept = delete;
};

//!\brief A bitvector representing the result of a call to `bulk_contains` of the seqan3::interleaved_bloom_filter.
template <typename FingerprintType>
class interleaved_xor_filter<FingerprintType>::membership_agent::binning_bitvector
{
private:
    //!\brief The underlying datatype to use.
    using data_type = sdsl::bit_vector;
    //!\brief The bitvector.
    data_type data{};

    friend class membership_agent;

    template <std::integral value_t>
    friend class counting_agent_type;

    template <std::integral value_t>
    friend class counting_vector;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    binning_bitvector() = default; //!< Defaulted.
    binning_bitvector(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector(binning_bitvector &&) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector &&) = default; //!< Defaulted.
    ~binning_bitvector() = default; //!< Defaulted.

    //!\brief Construct with given size.
    explicit binning_bitvector(size_t const size) :
        data(size)
    {}
    //!\}

    //!\brief Returns the number of elements.
    size_t size() const noexcept
    {
        return data.size();
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator to the first element of the container.
    auto begin() noexcept
    {
        return data.begin();
    }

    //!\copydoc begin()
    auto begin() const noexcept
    {
        return data.begin();
    }

    //!\brief Returns an iterator to the element following the last element of the container.
    auto end() noexcept
    {
        return data.end();
    }

    //!\copydoc end()
    auto end() const noexcept
    {
        return data.end();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Test for equality.
    friend bool operator==(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return lhs.data == rhs.data;
    }

    //!\brief Test for inequality.
    friend bool operator!=(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
     //!\brief Return the i-th element.
    auto operator[](size_t const i) noexcept
    {
        assert(i < size());
        return data[i];
    }

    //!\copydoc operator[]()
    auto operator[](size_t const i) const noexcept
    {
        assert(i < size());
        return data[i];
    }

    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}
};


/*!\brief Manages counting ranges of values for the seqan3::interleaved_xor_filter.
 * \attention Calling seqan3::interleaved_xor_filter::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_agent.cpp
 */
template <typename FingerprintType>
template <std::integral value_t>
class interleaved_xor_filter<FingerprintType>::counting_agent_type
{
private:
    //!\brief The type of the augmented seqan3::interleaved_xor_filter.
    using ixf_t = interleaved_xor_filter<FingerprintType>;

    //!\brief A pointer to the augmented seqan3::interleaved_xor_filter.
    ixf_t const * ixf_ptr{nullptr};

    //!\brief Store a seqan3::interleaved_xor_filter::membership_agent to call `bulk_contains`.
    typename ixf_t::membership_agent membership_agent;

public:
    class counting_vector;
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default; //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default; //!< Defaulted.
    ~counting_agent_type() = default; //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan3::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan3::interleaved_bloom_filter.
     */
    explicit counting_agent_type(ixf_t const & ixf) :
        ixf_ptr(std::addressof(ixf)), membership_agent(ixf), result_buffer(ixf.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    counting_vector result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::input_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector const & bulk_count(value_range_t && values) & noexcept
    {
        assert(ixf_ptr != nullptr);
        assert(result_buffer.size() == ixf_ptr->bin_count());

        static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, 0);

        for (auto && value : values)
            result_buffer += membership_agent.bulk_contains(value);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector const & bulk_count(value_range_t && values) && noexcept = delete;
    //!\}

};

/*!\brief A data structure that behaves like a std::vector and can be used to consolidate the results of multiple calls
 *        to seqan3::interleaved_xor_filter::membership_agent::bulk_contains.
 * \ingroup search_dream_index
 * \tparam value_t The type of the count. Must model std::integral.
 *
 * \details
 *
 * When using the seqan3::interleaved_xor_filter::membership_agent::bulk_contains operation, a common use case is to
 * add up, for example, the results for all k-mers in a query. This yields, for each bin, the number of k-mers of a
 * query that are in the respective bin. Such information can be used to apply further filtering or abundance estimation
 * based on the k-mer counts.
 *
 * The seqan3::counting_vector offers an easy way to add up the individual
 * seqan3::interleaved_xor_filter::membership_agent::binning_bitvector by offering an `+=` operator.
 *
 * The `value_t` template parameter should be chosen in a way that no overflow occurs if all calls to `bulk_contains`
 * return a hit for a specific bin. For example, `uint8_t` will suffice when processing short Illumina reads, whereas
 * long reads will require at least `uint32_t`.
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_vector.cpp
 */
template <typename FingerprintType>
template<std::integral value_t>
class interleaved_xor_filter<FingerprintType>::counting_agent_type<value_t>::counting_vector : public std::vector<value_t>
{
private:
    //!\brief The base type.
    using base_t = std::vector<value_t>;
    using ixf_t = interleaved_xor_filter<FingerprintType>;

    friend class counting_agent_type<value_t>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_vector() = default; //!< Defaulted.
    counting_vector(counting_vector const &) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector const &) = default; //!< Defaulted.
    counting_vector(counting_vector &&) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector &&) = default; //!< Defaulted.
    ~counting_vector() = default; //!< Defaulted.

    using base_t::base_t;
    typename ixf_t::membership_agent membership_agent;
    //!\}

    /*!\brief Bin-wise adds the bits of a seqan3::interleaved_xor_filter::membership_agent::binning_bitvector.
     * \tparam rhs_t The type of the right-hand side.
     *         Must be seqan3::interleaved_xor_filter::membership_agent::binning_bitvector.
     * \param rhs The seqan3::interleaved_xor_filter::membership_agent::binning_bitvector.
     * \attention The counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    template <typename rhs_t>
/*    //!\cond
        requires std::same_as<rhs_t, membership_agent::binning_bitvector>
    //!\endcond
*/
    counting_vector & operator+=(rhs_t const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        // Each iteration can handle 64 bits, so we need to iterate `((rhs.size() + 63) >> 6` many times
        for (size_t batch = 0, bin = 0; batch < ((rhs.size() + 63) >> 6); bin = 64 * ++batch)
        {
            size_t tmp = rhs.data.get_int(batch * 64); // get 64 bits starting at position `batch * 64`
            if (tmp ^ (1ULL<<63)) // This is a special case, because we would shift by 64 (UB) in the while loop.
            {
                while (tmp > 0)
                {
                    // Jump to the next 1 and increment the corresponding vector entry.
                    uint8_t step = std::countr_zero(tmp);
                    bin += step++;
                    tmp >>= step;
                    ++(*this)[bin++];
                }
            }
            else
            {
                ++(*this)[bin + 63];
            }
        }
        return *this;
    }

    /*!\brief Bin-wise addition of two `seqan3::counting_vector`s.
     * \param rhs The other seqan3::counting_vector.
     * \attention The seqan3::counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    counting_vector & operator+=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<value_t>());

        return *this;
    }
};

} // namespace seqan3
