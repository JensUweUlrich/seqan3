// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::scoring_scheme_concept.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

/*!\interface seqan3::scoring_scheme_concept <>
 * \brief A concept that requires that type be able to score two letters.
 * \tparam t            The type the concept check is performed on (the putative scoring scheme).
 * \tparam alphabet_t   The type of the first letter that you wish to score; must model seqan3::alphabet_concept.
 * \tparam alphabet2_t  The type of the second letter that you wish to score; must model seqan3::alphabet_concept;
 *                      defaults to `alphabet_t`.
 * \ingroup scoring
 *
 * \details
 *
 * This concept makes no assumptions about configurability or assignability of the scoring scheme, only the
 * ability to score the two letters is required.
 *
 */
/*!\name Requirements for seqan3::scoring_scheme_concept
 * \brief You can expect these members on all types that implement seqan3::scoring_scheme_concept.
 * \memberof seqan3::scoring_scheme_concept
 * \{
 */
/*!\typedef     typedef IMPLEMENTATION_DEFINED score_type;
 * \brief       The type returned by seqan3::scoring_scheme_concept::score(), usually a seqan3::arithmetic_concept.
 * \memberof seqan3::scoring_scheme_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual typedef (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn          score_type score(alph1_t const alph1, alph2_t const alph2);
 * \brief       Compute the score of two letters
 * \param alph1 First letter.
 * \param alph2 Second letter.
 * \memberof seqan3::scoring_scheme_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
//!\cond
template <typename t, alphabet_concept alphabet_t, alphabet_concept alphabet2_t = alphabet_t>
concept scoring_scheme_concept = requires (t scheme,
                                                alphabet_t const alph1,
                                                alphabet2_t const alph2)
{
    { scheme.score(alph1, alph2) } -> std::CommonReference<typename std::remove_reference_t<t>::score_type>;
    { scheme.score(alphabet_t{}, alphabet2_t{}) }
        -> std::CommonReference<typename std::remove_reference_t<t>::score_type>;
};
//!\endcond

} // namespace seqan3
