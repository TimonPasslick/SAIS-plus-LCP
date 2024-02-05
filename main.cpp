#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// I have a Windows laptop and the std library doesn't provide peak memory usage.
// #define LINUX_MEMORY_PEAK

#ifdef LINUX_MEMORY_PEAK
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#endif

template <typename L>
void print(const L& l)
{
	std::cout << std::setfill(' ');
	for (const auto i : l)
		if (i == -1)
			std::cout << "   ";
		else
			std::cout << std::setw(2) << i + 1 << ' ';
	std::cout << std::endl;
}

constexpr bool estimate_memory_peak {true};
long long memory_peak {0};
template <typename T, typename Int>
inline void allocated(const Int amount)
{
	if constexpr (estimate_memory_peak)
	{
		memory_peak += (long long) (sizeof(T) * amount);
	}
}

// SAIS algorithm
template <typename I, typename S>
std::vector<I> bucket_boundaries(const S& text, const I alphabet_size)
{
	std::vector<I> result(alphabet_size + 1, 0);
	allocated<I>(alphabet_size + 1); // n/2 in first recursion, so n total

	// histogram
	for (const auto character : text)
	{
		++result[character];
	}

	// prefix sum
	I sum {0};
	for (I& entry : result)
	{
		I copy {entry};
		entry = sum;
		sum += copy;
	}

	return result;
}
template <typename I, typename S>
void induce(std::vector<I>& suffix_array, std::vector<I>& inserted, const S& text, const std::vector<I>& bucket_bounds)
{
	std::fill(inserted.begin(), inserted.end(), 0);
	{ // scan left to right for L-type suffixes
		for (const I entry : suffix_array)
		{
			if (entry <= 0) // empty = -1 < 0
				continue;
			const auto candidate {text[entry - 1]};
			/* equivalent to if entry - 1 is L because:
			 *     if entry is S* then text[entry - 1] != text[entry]
			 *     entry can't be S, but not S* (not inserted yet)
			 */
			if (candidate >= text[entry])
				suffix_array[bucket_bounds[candidate]
				             + inserted[candidate]++] = entry - 1;
		}
	}
	std::fill(inserted.begin(), inserted.end(), 0);
	{ // scan right to left for S-type suffixes
		bool entry_is_s = true;
		auto bucket = (typename S::value_type)(bucket_bounds.size() - 2);
		const auto rend = suffix_array.rend();
		for (auto it = suffix_array.rbegin(); it != rend; ++it)
		{
			const auto one_based = rend - it;
			if (one_based == bucket_bounds[bucket])
			{
				--bucket;
				entry_is_s = true;
			}
			if (one_based == bucket_bounds[bucket + 1] - inserted[bucket])
				entry_is_s = false;
			if (*it <= 0) // empty = -1 < 0
				continue;
			const I candidate {text[*it - 1]};
			const auto right = text[*it];
			if (candidate < right || (entry_is_s && candidate == right))
				suffix_array[bucket_bounds[candidate + 1]
				             - ++inserted[candidate]] = *it - 1;
		}
	}
}
// i needs to be int16_t, int32_t or int64_t
template <typename I, typename S>
std::vector<I> get_suffix_array(const S& text, const I alphabet_size = 256)
{
	constexpr I empty {-1};

	std::vector<I> suffix_array(text.size(), empty);
	allocated<I>(text.size()); // n in base call, so 2n total
	const auto bucket_bounds = bucket_boundaries<I>(text, alphabet_size);
	std::vector<I> inserted(alphabet_size, 0);
	allocated<I>(alphabet_size); // n/2 in first recursion, so n total

	std::vector<I> lms;
	lms.reserve(text.size() / 2); // max amount
	allocated<I>(text.size() / 2); // n/2 in base call, so n total
	{ // put LMS-suffixes in text order at the end of buckets
		auto previous_character = *text.rbegin();
		bool previous_character_is_s {true};
		for (I i(text.size() - 1); i >= 0; --i)
		{
			const auto character = text[i];
			const bool is_s {
				character < previous_character
				|| (previous_character_is_s && character == previous_character)};
			if (previous_character_is_s && !is_s) // i + 1 is S*
			{
				lms.push_back(i + 1);

				const I bucket_end {bucket_bounds[previous_character + 1]};
				const I sa_pos = bucket_end - ++inserted[previous_character];
				suffix_array[sa_pos] = i + 1;
			}
			previous_character = character;
			previous_character_is_s = is_s;
		}
	}

	induce(suffix_array, inserted, text, bucket_bounds);

	// compute T': LMS substring ranks in text order
	// ranks are 1 based so that the sentinel can be 0
	std::vector<I> ranks(text.size(), 0); // gets shrinked to fit later on
	bool recursion_required {false};
		I rank_max {1};
	{ // insert ranks in text order with gaps
		const I text_end(text.size() - 1);
		ranks[text_end] = 1;
		I previous_text_index {text_end};
		for (I bucket {1}; bucket != alphabet_size; ++bucket)
		{
			const I end {bucket_bounds[bucket + 1]};
			bool in_lms_area {false};
			for (I i = end - inserted[bucket]; i != end; ++i)
			{
				const I text_start_index {suffix_array[i]};
				if (!in_lms_area && (text_start_index == 0 || text[text_start_index - 1] <= text[text_start_index]))
					continue; // normal S, not S*
				in_lms_area = true;
				I text_index {text_start_index};
				bool is_l {false};
				while (true)
				{ // find out if the LMS substrings are equal, if not, ++rank_max
					const auto current_character = text[text_index];
					if (text[previous_text_index] != current_character)
					{
						++rank_max;
						break;
					}
					const auto next_character = text[text_index + 1];
					if (is_l)
					{
						if (current_character < next_character)
						{ // arrived at last S*, LMS substrings equal
							recursion_required = true;
							break;
						}
					}
					else if (current_character > next_character)
						is_l = true;
					++previous_text_index;
					++text_index;
				}
				ranks[text_start_index] = rank_max;
				previous_text_index = text_start_index;
			}
		}
	}
	if (!recursion_required)
		return suffix_array;
	{ // remove gaps
		I size {0};
		for (const I rank : ranks)
			if (rank != 0)
				ranks[size++] = rank;
		ranks[size++] = 0; // sentinel
		ranks.resize(size);
		allocated<I>(size); // n/2 in base call, so n total
	}
	const auto order = get_suffix_array(ranks, rank_max + 1);
	std::fill(suffix_array.begin(), suffix_array.end(), empty);
	std::fill(inserted.begin(), inserted.end(), 0);
	for (auto it {order.crbegin()}; it < order.crend() - 1; ++it)
	{
		const I text_index {*(lms.rbegin() + *it)};
		const I character {text[text_index]};
		suffix_array[bucket_bounds[character + 1]
		             - ++inserted[character]] = text_index;
	}
	induce(suffix_array, inserted, text, bucket_bounds);
	return suffix_array;
}

template <typename I>
const std::vector<I> get_lcp_array_naive(const std::vector<std::uint8_t>& text, const std::vector<I>& suffix_array)
{
	std::vector<I> lcp_array;
	lcp_array.reserve(suffix_array.size());
	lcp_array.push_back(0);
	for (auto it = suffix_array.begin() + 1; it != suffix_array.end(); ++it)
	{
		const I i_start {*(it - 1)};
		I i {i_start};
		I j {*it};
		while (i != text.size() && j != text.size() && text[i] == text[j])
		{
			++i;
			++j;
		}
		lcp_array.push_back(i - i_start);
	}
	return lcp_array;
}

template <typename I>
const std::vector<I> get_isa(const std::vector<I>& sa)
{
	std::vector<I> isa(sa.size());
	for (I i {0}; i != sa.size(); ++i)
		isa[sa[i]] = i;
	return isa;
}
template <typename I>
const std::vector<I> get_lcp_array_kasai(const std::vector<std::uint8_t>& t, const std::vector<I>& sa)
{
	const auto isa = get_isa(sa);
	std::vector<I> lcp(sa.size(), 0);
	I L = 0;
	for (int i {0}; i < sa.size(); ++i)
	{
		if (isa[i] == 0)
			break;
		I j {sa[isa[i] - 1]};
		while (t[i + L] == t[j + L])
			++L;
		lcp[isa[i]] = L;
		--L;
		if (L < 0)
			L = 0;
	}
	return lcp;
}

template <typename I>
const std::vector<I> get_lcp_array_phi(const std::vector<std::uint8_t>& t, const std::vector<I>& sa)
{
	std::vector<I> phi(sa.size());
	const I n = sa.size() - 1;
	phi[n] = sa[n];
	for (I i {1}; i != sa.size(); ++i)
		phi[sa[i]] = sa[i - 1];
	I L {0};
	for (I i{0}; i != sa.size(); ++i)
	{
		I j {phi[i]};
		while (t[i + L] == t[j + L])
			++L;
		phi[i] = L;
		--L;
		if (L < 0)
			L = 0;
	}
	std::vector<I> lcp;
	lcp.reserve(sa.size());
	for (I i{0}; i != sa.size(); ++i)
		lcp.push_back(phi[sa[i]]);
	return sa;
}

// copied and adapted Jan 27th from:
// https://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
std::vector<std::uint8_t> get_file_contents(const char *filename)
{
	std::ifstream in {filename, std::ios::in | std::ios::binary};
	if (!in)
	{
		std::cerr << "File could not be opened." << std::endl;
		std::exit(1);
	}
	in.seekg(0, std::ios::end);
	std::size_t length {static_cast<std::size_t>(in.tellg())};
	std::vector<std::uint8_t> contents;
	contents.reserve(length + 1); // +1 for sentinel
	contents.resize(length);
	in.seekg(0, std::ios::beg);
	in.read(reinterpret_cast<char*>(&contents[0]), contents.size());
	in.close();
	if (std::find(contents.begin(), contents.end(), 0) != contents.end())
	{
		std::cerr << "Null characters not supported." << std::endl;
		std::exit(1);
	}
	contents.push_back(0); // sentinel required for SA and LCP array construction
	return contents;
}

std::int64_t get_execution_time(const std::function<void()>& function)
{
	using namespace std::chrono;
	const auto before = high_resolution_clock::now();
	function();
	return duration_cast<milliseconds>(high_resolution_clock::now() - before).count();
}

#ifdef LINUX_MEMORY_PEAK
long long get_memory_peak()
{
	char* buffer;
	std::size_t size;
	FILE* file {open_memstream(&buffer, &size)};
	malloc_info(0, file);
	const std::string info {buffer};
	free(buffer);
	const std::string prefix {"<system type=\"max\" size=\""};
	const std::size_t prefix_position {info.find(prefix)};
	const std::size_t number_start {prefix_position + prefix.size()};
	char* number_end;
	return std::strtoll(info.data() + number_start, &number_end, 10);
}
#endif

template <typename I>
void run(std::int64_t& sa_construction_time,
         std::int64_t& lcp_naive_construction_time,
         std::int64_t& lcp_kasai_construction_time,
         std::int64_t& lcp_phi_construction_time,
		 const std::vector<std::uint8_t>& text)
{
	std::vector<I> suffix_array;
	sa_construction_time = get_execution_time([&]() {
		suffix_array = get_suffix_array<I>(text);
	});
	lcp_naive_construction_time = get_execution_time([&]() {
		get_lcp_array_naive<I>(text, suffix_array);
	});
	lcp_kasai_construction_time = get_execution_time([&]() {
		get_lcp_array_kasai<I>(text, suffix_array);
	});
	lcp_phi_construction_time = get_execution_time([&]() {
		get_lcp_array_phi<I>(text, suffix_array);
	});
}
int main(int argc, char** argv)
{
	std::vector<std::uint8_t> text;
	switch (argc)
	{
	case 1:
		text = {'m', 'i', 's', 's', 'i', 's', 's', 'i', 'p', 'p', 'i', 0};
		break;
	case 2:
		text = get_file_contents(argv[1]);
		break;
	default:
		std::cerr << "too many arguments, only input file expected" << std::endl;
		return 1;
	}
	std::int64_t sa_construction_time;
	std::int64_t lcp_naive_construction_time;
	std::int64_t lcp_kasai_construction_time;
	std::int64_t lcp_phi_construction_time;
	if (text.size() <= std::numeric_limits<std::int16_t>::max())
	{
		run<std::int16_t>(
			sa_construction_time,
			lcp_naive_construction_time,
			lcp_kasai_construction_time,
			lcp_phi_construction_time,
			text
		);
	}
	else if (text.size() <= std::numeric_limits<std::int32_t>::max())
	{
		run<std::int32_t>(
			sa_construction_time,
			lcp_naive_construction_time,
			lcp_kasai_construction_time,
			lcp_phi_construction_time,
			text
		);
	}
	else if (text.size() <= std::numeric_limits<std::int64_t>::max())
	{
		run<std::int64_t>(
			sa_construction_time,
			lcp_naive_construction_time,
			lcp_kasai_construction_time,
			lcp_phi_construction_time,
			text
		);
	}
#ifdef LINUX_MEMORY_PEAK
	memory_peak = get_memory_peak();
#endif
	// Print the factor to see how much larger the memory needed is than the text.
	// Worst case for large texts: 6 * sizeof(std::intxx_t)
	const auto factor = (double) memory_peak / text.size();
	constexpr long long megabyte {2 << 20};
	// rounded division
	memory_peak = (memory_peak + megabyte / 2) / megabyte;
	std::cout << "RESULT name=TimonPasslick"
		<< " sa_construction_time=" << sa_construction_time
		<< " sa_construction_memory=" << memory_peak
		<< " lcp_naive_construction_time=" << lcp_naive_construction_time
		<< " lcp_kasai_construction_time=" << lcp_kasai_construction_time
		<< " lcp_phi_construction_time=" << lcp_phi_construction_time
		<< std::endl;
}
