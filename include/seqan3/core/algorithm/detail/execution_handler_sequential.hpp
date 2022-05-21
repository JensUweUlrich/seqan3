// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::execution_handler_sequential.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <seqan3/std/ranges>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Handles the sequential execution of algorithms.
 * \ingroup core_algorithm
 *
 * \details
 *
 * This execution handler implements a blocking execute. This means a call to
 * seqan3::detail::execution_handler_sequential::execute blocks until the call to the algorithm finished.
 * This handler can be used in combination with the seqan3::detail::algorithm_executor_blocking to invoke the
 * algorithms on the given algorithm input.
 */
struct execution_handler_sequential
{
public:

    /*!\brief Executes the algorithm with the given input and callback.
     * \tparam algorithm_t The type of the algorithm.
     * \tparam algorithm_input_t The input type to invoke the algorithm with.
     * \tparam callback_t The type of the callable invoked by the algorithm after generating a new result.
     *
     * \param[in] algorithm The algorithm to invoke.
     * \param[in] input The input of the algorithm.
     * \param[in] callback A callable which will be invoked on each result generated by the algorithm.
     */
    template <typename algorithm_t, typename algorithm_input_t, typename callback_t>
        requires std::invocable<algorithm_t, algorithm_input_t, callback_t>
    void execute(algorithm_t && algorithm, algorithm_input_t && input, callback_t && callback)
    {
        algorithm(std::forward<algorithm_input_t>(input), std::forward<callback_t>(callback));
    }

    /*!\brief Sequentially executes the algorithm for every element of the given input range.
     * \tparam algorithm_t The type of the algorithm.
     * \tparam algorithm_input_range_t The input range type.
     * \tparam callback_t The type of the callable invoked by the algorithm after generating a new result.
     *
     * \param[in] algorithm The algorithm to invoke.
     * \param[in] input_range The input range to process sequentially.
     * \param[in] callback A callable which will be invoked on each result generated by the algorithm for a given input.
     *
     * \details
     *
     * Effectively calls seqan3::detail::execution_handler_sequential::execute on every element of the given input
     * range.
     */
    template <std::copy_constructible algorithm_t,
              std::ranges::input_range algorithm_input_range_t,
              std::copy_constructible callback_t>
        requires std::invocable<algorithm_t, std::ranges::range_reference_t<algorithm_input_range_t>, callback_t>
    void bulk_execute(algorithm_t && algorithm, algorithm_input_range_t && input_range, callback_t && callback)
    {
        for (auto && input : input_range)
            execute(algorithm, std::forward<decltype(input)>(input), callback);

        wait();
    }

    //!\brief Waits for the submitted jobs to finish (noop).
    void wait() noexcept
    {
        // noop
    }
};

} // namespace seqan3
