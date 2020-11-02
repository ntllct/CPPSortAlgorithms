#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "sort.hpp"

#ifdef __linux__
#include <cxxabi.h>
#elif _WIN32
#include <windows.h>
#else
// reserved
#endif
// for random string and shuffle
std::mt19937_64 random64{ 1234 };
// Default array size
size_t array_size = 50'000;
// Maximum random string length
const size_t MAX_STRING_LENGTH = 15;  // + 1

enum class shuffletype : unsigned char {
    random,
    sorted,
    swap_first_and_last,
    reversed,
    almost_sorted_4,
    almost_sorted_8,
    almost_sorted_16,
};
std::vector<std::pair<std::string, shuffletype>> shufflenames = {
    {"Sort random values", shuffletype::random},
    {"Sort sorted values", shuffletype::sorted},
    {"Sort almost sorted values", shuffletype::swap_first_and_last},
    {"Sort reversed sorted values", shuffletype::reversed},
    {"Sort almost sorted values  4", shuffletype::almost_sorted_4},
    {"Sort almost sorted values  8", shuffletype::almost_sorted_8},
    {"Sort almost sorted values 16", shuffletype::almost_sorted_16} };
// Shuffle's names in command line arguments
std::map<std::string, shuffletype> shuffleames_in_cmd = {
    {"random", shuffletype::random},
    {"sorted", shuffletype::sorted},
    {"flswap", shuffletype::swap_first_and_last},
    {"reversed", shuffletype::reversed},
    {"s4", shuffletype::almost_sorted_4},
    {"s8", shuffletype::almost_sorted_8},
    {"s16", shuffletype::almost_sorted_16} };

enum class sortalgorithm : unsigned char {
    bubble,
    oddeven,
    shaker,
    comb,
    insertion,
    selection,
    shell,
    shellhib,
    shellpratt,
    shellsedgwick,
    shell1,
    shell2,
    shell3,
    tree,
    gnome,
    heap,
    quick,
    merge,
    bucket,
    lsd,
    msd,
    bitonic,
    tim,
};
std::vector<std::tuple<std::string, sortalgorithm, bool>> sortnames = {
    {"Bubble sort", sortalgorithm::bubble, true},
    {"Odd-even sort", sortalgorithm::oddeven, true},
    {"Shaker sort", sortalgorithm::shaker, true},
    // false -> skip big array because it is too slow
    {"Comb sort", sortalgorithm::comb, false},
    {"Insertion sort", sortalgorithm::insertion, true},
    {"Selection sort", sortalgorithm::selection, true},
    {"Shell sort", sortalgorithm::shell, false},
    {"Shellhib sort", sortalgorithm::shellhib, false},
    {"Shellpratt sort", sortalgorithm::shellpratt, false},
    {"Shellsedgwick sort", sortalgorithm::shellsedgwick, false},
    {"Shell1 sort", sortalgorithm::shell1, false},
    {"Shell2 sort", sortalgorithm::shell2, false},
    {"Shell3 sort", sortalgorithm::shell3, false},
    {"Tree sort", sortalgorithm::tree, false},
    {"Gnome sort", sortalgorithm::gnome, true},
    {"Heap sort", sortalgorithm::heap, false},
    {"Quick sort", sortalgorithm::quick, false},
    {"Merge sort", sortalgorithm::merge, false},
    {"Bucket sort", sortalgorithm::bucket, false},
    {"Lsd sort", sortalgorithm::lsd, false},
    {"Msd sort", sortalgorithm::msd, false},
    {"Bitonic sort", sortalgorithm::bitonic, false},
    {"Tim sort", sortalgorithm::tim, false} };
// For recognizing from console command
std::map<std::string, std::pair<std::string, sortalgorithm>>
sortnames_in_cmd = {
    {"bubble", {"Bubble sort", sortalgorithm::bubble}},
    {"oddeven", {"Odd-even sort", sortalgorithm::oddeven}},
    {"shaker", {"Shaker sort", sortalgorithm::shaker}},
    {"comb", {"Comb sort", sortalgorithm::comb}},
    {"insertion", {"Insertion sort", sortalgorithm::insertion}},
    {"selection", {"Selection sort", sortalgorithm::selection}},
    {"shell", {"Shell sort", sortalgorithm::shell}},
    {"shellhib", {"Shellhib sort", sortalgorithm::shellhib}},
    {"shellpratt", {"Shellpratt sort", sortalgorithm::shellpratt}},
    {"shellsedwick", {"Shellsedgwick sort", sortalgorithm::shellsedgwick}},
    {"shell1", {"Shell1 sort", sortalgorithm::shell1}},
    {"shell2", {"Shell2 sort", sortalgorithm::shell2}},
    {"shell3", {"Shell3 sort", sortalgorithm::shell3}},
    {"tree", {"Tree sort", sortalgorithm::tree}},
    {"gnome", {"Gnome sort", sortalgorithm::gnome}},
    {"heap", {"Heap sort", sortalgorithm::heap}},
    {"quick", {"Quick sort", sortalgorithm::quick}},
    {"merge", {"Merge sort", sortalgorithm::merge}},
    {"bucket", {"Bucket sort", sortalgorithm::bucket}},
    {"lsd", {"Lsd sort", sortalgorithm::lsd}},
    {"msd", {"Msd sort", sortalgorithm::msd}},
    {"bitonic", {"Bitonic sort", sortalgorithm::bitonic}},
    {"tim", {"Tim sort", sortalgorithm::tim}} };

template <typename T>
std::string type_name() {
    int status{ 0 };
    std::string tname = typeid(T).name();
#ifdef __linux__
    char* demangled_name =
        abi::__cxa_demangle(tname.c_str(), nullptr, nullptr, &status);
    if (status == 0) {
        tname = demangled_name;
        std::free(demangled_name);
    }
#elif _WIN32
    // reserved
#else
    // reserved
#endif
    return (tname);
}

template <typename T>
void shuffle_data(std::vector<T>& data, shuffletype st) {
    switch (st) {
    case shuffletype::sorted:
        break;
    case shuffletype::swap_first_and_last:
        if (data.size() > 1) std::swap(data.front(), data.back());
        break;
    case shuffletype::reversed:
        std::reverse(data.begin(), data.end());
        break;
    case shuffletype::almost_sorted_4:
        for (auto it = data.begin(); it < data.end(); it += 4) {
            auto end = (it + 4 < data.end()) ? it + 4 : data.end();
            std::shuffle(it, end, random64);
        }
        break;
    case shuffletype::almost_sorted_8:
        for (auto it = data.begin(); it < data.end(); it += 8) {
            auto end = (it + 8 < data.end()) ? it + 8 : data.end();
            std::shuffle(it, end, random64);
        }
        break;
    case shuffletype::almost_sorted_16:
        for (auto it = data.begin(); it < data.end(); it += 16) {
            auto end = (it + 16 < data.end()) ? it + 16 : data.end();
            std::shuffle(it, end, random64);
        }
        break;
    default:
        std::shuffle(data.begin(), data.end(), random64);
        break;
    }
}

template <typename T>
std::vector<T> gen(size_t length, size_t nUniques, shuffletype st) {
    std::vector<T> out(length, 0);
    if (length == 0) return (out);
    size_t uniques = nUniques == 0 ? 1 : length / nUniques;
    if (uniques == 0) uniques++;

    if constexpr (std::is_floating_point<T>::value) {
        for (size_t i = 0; i < length; ++i) {
            out[i] = nUniques == 0 ? static_cast<T>(i) : static_cast<T>(i / uniques);
            out[i] -= static_cast<T>(length / 2);
        }
    }
    else if constexpr (std::is_signed<T>::value) {
        for (size_t i = 0; i < length; ++i) {
            out[i] = nUniques == 0 ? static_cast<T>(i) : static_cast<T>(i / uniques);
            out[i] -= static_cast<T>(length / 2);
        }
    }
    else if constexpr (std::is_unsigned<T>::value) {
        for (size_t i = 0; i < length; ++i) {
            out[i] = nUniques == 0 ? static_cast<T>(i) : static_cast<T>(i / uniques);
        }
    }
    std::sort(out.begin(), out.end());
    shuffle_data(out, st);
    return (out);
}
template <>
std::vector<std::string> gen<std::string>(size_t length, size_t nUniques,
    shuffletype st) {
    std::vector<std::string> out(length, "");
    if (length == 0) return (out);

    auto random_string = []() -> std::string {
        static const char letters[] = "abcdefghijklmnopqrstuvwxyz";
        size_t nLetters = random64() % MAX_STRING_LENGTH + 1;
        std::string result(nLetters, 'x');
        for (size_t i = 0; i < nLetters; ++i)
            result[i] = letters[random64() % (sizeof(letters) - 1)];
        return (result);
    };

    size_t last_unique = 0;
    out[0] = random_string();
    size_t uniques = nUniques == 0 ? 1 : length / nUniques;
    if (uniques == 0) uniques++;
    for (size_t i = 1; i < length; ++i) {
        if (nUniques == 0) {
            out[i] = random_string();
        }
        else {
            if (i / uniques != last_unique) {
                last_unique = i / uniques;
                out[i] = random_string();
            }
            else {
                out[i] = out[i - 1];
            }
        }
    }
    std::sort(out.begin(), out.end());
    shuffle_data(out, st);
    return (out);
}

// Command line parser
class CLParser {
public:
    CLParser(int argc, const char* const* argv) {
        std::copy(argv + 1, argv + argc, std::back_inserter(tokens));
    }
    template <typename T>
    std::string get(T&& option) const noexcept {
        auto it = std::find(tokens.begin(), tokens.end(), std::forward<T>(option));
        if (it != tokens.end() && ++it != tokens.end()) return (*it);
        return (std::string());
    }
    template <typename T>
    bool check(T&& option) const noexcept {
        auto it = std::find(tokens.begin(), tokens.end(), std::forward<T>(option));
        return (it != tokens.end());
    }

private:
    std::vector<std::string> tokens{};
};

// Default array sizes
const std::vector<size_t> nValues{ 500000, 50000, 5000, 500, 50, 5, 2, 1, 0 };

template <typename T>
bool sort_case(std::vector<T>& data, sortalgorithm algo) {
    switch (algo) {
    case sortalgorithm::bubble:
        sort::bubble(data.begin(), data.end());
        break;
    case sortalgorithm::oddeven:
        sort::oddeven(data.begin(), data.end());
        break;
    case sortalgorithm::shaker:
        sort::shaker(data.begin(), data.end());
        break;
    case sortalgorithm::comb:
        sort::comb(data.begin(), data.end());
        break;
    case sortalgorithm::insertion:
        sort::insertion(data.begin(), data.end());
        break;
    case sortalgorithm::selection:
        sort::selection(data.begin(), data.end());
        break;
    case sortalgorithm::shell:
        sort::shell(data.begin(), data.end());
        break;
    case sortalgorithm::shellhib:
        sort::shellhib(data.begin(), data.end());
        break;
    case sortalgorithm::shellpratt:
        sort::shellpratt(data.begin(), data.end());
        break;
    case sortalgorithm::shellsedgwick:
        sort::shellsedgwick(data.begin(), data.end());
        break;
    case sortalgorithm::shell1:
        sort::shell1(data.begin(), data.end());
        break;
    case sortalgorithm::shell2:
        sort::shell2(data.begin(), data.end());
        break;
    case sortalgorithm::shell3:
        sort::shell3(data.begin(), data.end());
        break;
    case sortalgorithm::tree:
        sort::tree(data.begin(), data.end());
        break;
    case sortalgorithm::gnome:
        sort::gnome(data.begin(), data.end());
        break;
    case sortalgorithm::quick:
        sort::quick(data.begin(), data.end());
        break;
    case sortalgorithm::heap:
        sort::heap(data.begin(), data.end());
        break;
    case sortalgorithm::merge:
        sort::merge(data.begin(), data.end());
        break;
    case sortalgorithm::bucket:
        if constexpr (!std::is_same<T, std::string>::value) {
            sort::bucket(data.begin(), data.end());
            break;
        }
        return (false);
    case sortalgorithm::lsd:
        if constexpr (std::is_unsigned<T>::value) {
            sort::lsd(data.begin(), data.end());
            break;
        }
        return (false);
    case sortalgorithm::msd:
        if constexpr (std::is_unsigned<T>::value) {
            sort::msd(data.begin(), data.end());
            break;
        }
        return (false);
    case sortalgorithm::bitonic:
        if constexpr (std::is_unsigned<T>::value) {
            sort::bitonic(data.begin(), data.end());
            break;
            break;
        }
        return (false);
    case sortalgorithm::tim:
        sort::tim(data.begin(), data.end());
        break;
    default:
        return (false);
        break;
    }
    return (true);
}

template <typename T>
void benchmark_sort(
    const std::tuple<std::string, sortalgorithm, bool>& sort_algo,
    shuffletype st) {
    std::vector<T> data{};
    std::cout << std::setw(20) << std::left << std::get<std::string>(sort_algo) << ": ";
    for (size_t i = 0; i < nValues.size(); ++i) {
        // Skip taugh tasks
        if (i == 0 && std::get<bool>(sort_algo)) {
            std::cout << std::right << std::setw(15) << "[disabled]" << std::flush;
            continue;
        }
        data = gen<T>(nValues[i], 0, st);
        auto start_point = std::chrono::system_clock::now();
        bool is_processed = sort_case<T>(data, std::get<sortalgorithm>(sort_algo));
        if (is_processed) {
            auto end_point = std::chrono::system_clock::now();
            std::chrono::duration<double> dur = end_point - start_point;
#ifdef __linux__
            if (std::is_sorted(data.begin(), data.end()))
                std::cout << "\033[32m";
            else
                std::cout << "\033[31m";
            std::cout << std::right << std::setw(15) << dur.count();
            std::cout << "\033[0m" << std::flush;
#elif _WIN32
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            if (std::is_sorted(data.begin(), data.end()))
                SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
            else
                SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
            std::cout << std::right << std::setw(15) << dur.count();
            SetConsoleTextAttribute(hConsole, 7);
            std::cout << std::flush;
#else
            // reserved
#endif
        }
        else {
            std::cout << std::right << std::setw(15) << "[disabled]" << std::flush;
        }
    }
    std::cout << std::endl;
}
// Test all sort algorithms
template <typename T>
void test_sort_algorithms() {
    for (const auto& t : shufflenames) {
        std::cout << std::endl
            << t.first << " (" << type_name<T>() << ")" << std::endl
            << std::endl;
        std::cout << std::setw(22) << std::left << "Array length: ";
        for (const auto n : nValues)
            std::cout << std::right << std::setw(15) << n;
        std::cout << std::endl << std::endl;
        for (const auto& a : sortnames) benchmark_sort<T>(a, t.second);
    }
}
// Test sort algorithm specified from console
template <typename T>
void test_sort_specified_algorithm(
    const std::pair<std::string, sortalgorithm>& sa, const std::string& st,
    const std::string& vt) {
    std::vector<T> data;
    data = gen<T>(array_size, 0, shuffleames_in_cmd[st]);

    std::cout << std::left << std::setw(15) << sa.first << '|';
    std::cout << std::left << std::setw(10) << st << '|';
    std::cout << std::left << std::setw(8) << vt << '|';
    std::cout << std::right << std::setw(8) << data.size() << ": ";
    std::cout << std::left;

    auto start_point = std::chrono::system_clock::now();
    bool is_processed = sort_case<T>(data, sa.second);
    if (is_processed) {
        auto end_point = std::chrono::system_clock::now();
        std::chrono::duration<double> dur = end_point - start_point;
        if (std::is_sorted(data.begin(), data.end()))
            std::cout << dur.count() << std::endl;
        else
            std::cout << "[error]" << std::endl;
    }
    else {
        std::cout << "[skip]" << std::endl;
    }
}
std::map<std::string, void (*)(const std::pair<std::string, sortalgorithm>&,
    const std::string&, const std::string&)>
    instanciated_functions{
        {"s8",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<int8_t>(sa, st, vt);
         }},
        {"u8",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<uint8_t>(sa, st, vt);
         }},
        {"s16",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<int16_t>(sa, st, vt);
         }},
        {"u16",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<uint16_t>(sa, st, vt);
         }},
        {"s32",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<int32_t>(sa, st, vt);
         }},
        {"u32",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<uint32_t>(sa, st, vt);
         }},
        {"s64",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<int64_t>(sa, st, vt);
         }},
        {"u64",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<uint64_t>(sa, st, vt);
         }},
        {"f",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<float>(sa, st, vt);
         }},
        {"d",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<double>(sa, st, vt);
         }},
        {"ld",
         [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<long double>(sa, st, vt);
         }},
        {"string", [](const auto& sa, const auto& st, const auto& vt) {
           test_sort_specified_algorithm<std::string>(sa, st, vt);
         }} };
void sort_assert(bool result, void (*fn)(const char*), const char* str) {
    if (result) return;
    fn(str);
    std::exit(EXIT_SUCCESS);
}

namespace console_messages {
    void print_help_and_exit(const char* str) {
        std::cout << "Usage: " << str << " [OPTIONS]" << std::endl;
        std::cout << "\t-h, --help                     print help" << std::endl;
        std::cout
            << "\t-a, --algo sortnames           comma separated algorithm's names"
            << std::endl;
        std::cout << "\t                               "
            "bubble,oddeven,shaker,comb,insertion,selection,shell,"
            << std::endl;
        std::cout << "\t                               "
            "shellhib,shellpratt,shellsedgwick,shell1,shell2,shell3,"
            << std::endl;
        std::cout << "\t                               "
            "tree,gnome,heap,quick,merge,bucket,lsd,msd,bitonic,tim"
            << std::endl;
        std::cout << "\t-s, --size N                   array size" << std::endl;
        std::cout << "\t-t, --type value types         comma separated type names"
            << std::endl;
        std::cout
            << "\t                               s8  - int8_t;      u8  - uint8_t"
            << std::endl;
        std::cout
            << "\t                               s16 - int16_t;     u16 - uint16_t"
            << std::endl;
        std::cout
            << "\t                               s32 - int32_t;     u32 - uint32_t"
            << std::endl;
        std::cout << "\t                               s64 - int64_t;     u64 - uint64_t "
            << std::endl;
        std::cout << "\t                               f - float; d - double; ld - "
            "long double"
            << std::endl;
        std::cout << "\t                               string - std::string"
            << std::endl;
        std::cout << "\t-sh, --shuffle shuffle_type    comma separated shuffle "
            "algorithm's types"
            << std::endl;
        std::cout << "\t                               "
            "sorted,reversed,random,flswap,s4,s8,s16"
            << std::endl;
        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "\t" << str
            << " -a comb,shell -s 10000 -t s8,u64 -sh sorted,reversed,random"
            << std::endl;
        std::cout << "\t" << str
            << " -a bubble,oddeven -s 10000 -t u8,ld -sh sorted,reversed,random"
            << std::endl;
        std::exit(EXIT_SUCCESS);
    }
    void algos_not_found(const char* str) {
        std::cout << "Error: sorting algorithms is not specified!" << std::endl;
        std::cout << "Try " << str << " -a bubble -s 10000 -t s8,u64 -sh sorted"
            << std::endl;
        std::cout << "    " << str << " -a bubble,comb -s 10000 -t s8,u64 -sh sorted"
            << std::endl;
        std::cout << "    " << str << " -h # for more information" << std::endl;
    }
    void array_size_is_not_specified(const char* str) {
        std::cout << "Error: array size is not specified!" << std::endl;
        std::cout << "Try " << str << " -s 10000" << std::endl;
        std::cout << "    " << str << " --size 5" << std::endl;
    }
    void array_size_is_not_valid(const char* str) {
        std::cout << "Error: array size is not valid!" << std::endl;
        std::cout << "Try " << str << " -s 100000" << std::endl;
        std::cout << "    " << str << " --size 5" << std::endl;
    }
    void variable_types_is_not_specified(const char* str) {
        std::cout << "Error: variable types is not specified!" << std::endl;
        std::cout << "Try " << str
            << " -t s32 -a bubble,comb -s 10000 -sh sorted,reversed,random"
            << std::endl;
        std::cout << "    " << str
            << " -t s32,u32 -a bubble,comb -s 10000 -sh sorted,reversed,random"
            << std::endl;
        std::cout << "    " << str << " -h # for more information" << std::endl;
    }
    void shuffle_types_is_not_specified(const char* str) {
        std::cout << "Error: shuffle types is not specified!" << std::endl;
        std::cout << "Try " << str << " -sh random -t s32,u32 -a bubble,comb -s 10000"
            << std::endl;
        std::cout << "    " << str
            << " -sh reversed,random -t s32,u32 -a bubble,comb -s 10000"
            << std::endl;
        std::cout << "    " << str << " -h # for more information" << std::endl;
    }
    void unknown_algorithm_type(const char* str) {
        std::cout << "Unknown algorithm type: " << str << std::endl;
    }
    void unknown_variable_type(const char* str) {
        std::cout << "Unknown variable type: " << str << std::endl;
    }
    void unknown_shuffle_type(const char* str) {
        std::cout << "Unknown shuffle type: " << str << std::endl;
    }
    void specify_3_options(const char* str) {
        std::cout << "Specify 3 options '-a', '-t' and '-sh' or do not specify any "
            "options at all."
            << std::endl;
        std::cout << "Example: " << str << " -a bubble,comb -s 10000 -t u16 -sh random"
            << std::endl;
    }
}  // namespace console_messages
namespace param_extractors {
    std::vector<std::pair<std::string, sortalgorithm>> get_algorithms(
        CLParser& clp, const char* programm_name) {
        std::vector<std::pair<std::string, sortalgorithm>> sortalgorithms{};
        if (clp.check("-a") || clp.check("--algo")) {
            std::string algos = clp.check("-a") ? clp.get("-a") : clp.get("--algo");
            sort_assert(!algos.empty(), console_messages::algos_not_found,
                programm_name);
            std::stringstream ss_algos(algos);
            while (ss_algos.good()) {
                std::string algo{};
                std::getline(ss_algos, algo, ',');
                if (algo.empty()) continue;
                auto it = sortnames_in_cmd.find(algo);
                sort_assert(it != sortnames_in_cmd.end(),
                    console_messages::unknown_algorithm_type, algo.c_str());
                sortalgorithms.push_back(it->second);
            }
        }
        return (sortalgorithms);
    }
    bool check_array_size(CLParser& clp, const char* programm_name) {
        if (clp.check("-s") || clp.check("--size")) {
            std::string size = clp.check("-s") ? clp.get("-s") : clp.get("--size");
            sort_assert(!size.empty(), console_messages::array_size_is_not_specified,
                programm_name);
            sort_assert(isdigit(size[0]) && size[0] != '0',
                console_messages::array_size_is_not_valid, programm_name);
            array_size = std::stoul(size);
            return (true);
        }
        return (false);
    }
    std::vector<std::string> get_variable_types(CLParser& clp,
        const char* programm_name) {
        std::vector<std::string> variable_types{};
        if (clp.check("-t") || clp.check("--type")) {
            std::string types = clp.check("-t") ? clp.get("-t") : clp.get("--type");
            sort_assert(!types.empty(),
                console_messages::variable_types_is_not_specified,
                programm_name);
            std::stringstream ss_types(types);
            while (ss_types.good()) {
                std::string type;
                std::getline(ss_types, type, ',');
                if (type.empty()) continue;
                sort_assert(instanciated_functions.count(type) > 0,
                    console_messages::unknown_variable_type, type.c_str());
                variable_types.emplace_back(std::move(type));
            }
        }
        return (variable_types);
    }
    std::vector<std::string> get_shuffle_types(CLParser& clp,
        const char* programm_name) {
        std::vector<std::string> shuffle_types{};
        if (clp.check("-sh") || clp.check("--shuffle")) {
            std::string types =
                clp.check("-sh") ? clp.get("-sh") : clp.get("--shuffle");
            sort_assert(!types.empty(),
                console_messages::shuffle_types_is_not_specified,
                programm_name);
            std::stringstream ss_types(types);
            while (ss_types.good()) {
                std::string type;
                std::getline(ss_types, type, ',');
                if (type.empty()) continue;
                sort_assert(shuffleames_in_cmd.count(type) > 0,
                    console_messages::unknown_shuffle_type, type.c_str());
                shuffle_types.emplace_back(std::move(type));
            }
        }
        return (shuffle_types);
    }
}  // namespace param_extractors

int main(int argc, const char* argv[]) {
    CLParser clp(argc, argv);
    if (clp.check("-h") || clp.check("--help"))
        console_messages::print_help_and_exit(argv[0]);
    // Load parameters from command line arguments
    auto sortalgorithms = param_extractors::get_algorithms(clp, argv[0]);
    auto is_size_specified = param_extractors::check_array_size(clp, argv[0]);
    auto variable_types = param_extractors::get_variable_types(clp, argv[0]);
    auto shuffle_types = param_extractors::get_shuffle_types(clp, argv[0]);
    // Setup console
    std::cout.precision(8);
    std::cout << std::fixed;
    // Some checks
    bool has_empty = sortalgorithms.empty() || variable_types.empty() || shuffle_types.empty();
    bool has_notempty = !sortalgorithms.empty() || !variable_types.empty() || !shuffle_types.empty();
    if (has_empty == has_notempty ||
        (is_size_specified && has_empty)) {
        console_messages::specify_3_options(argv[0]);
        std::exit(EXIT_FAILURE);
    }
    // Perform algorithms
    if (!sortalgorithms.empty() && !variable_types.empty() &&
        !shuffle_types.empty()) {
        for (const auto& sa : sortalgorithms)
            for (const auto& st : shuffle_types)
                for (const auto& vt : variable_types)
                    instanciated_functions[vt](sa, st, vt);
    }
    else {
        test_sort_algorithms<uint64_t>();
        test_sort_algorithms<uint32_t>();
        test_sort_algorithms<uint16_t>();
        test_sort_algorithms<uint8_t>();
        test_sort_algorithms<int64_t>();
        test_sort_algorithms<int32_t>();
        test_sort_algorithms<int16_t>();
        test_sort_algorithms<int8_t>();
        test_sort_algorithms<float>();
        test_sort_algorithms<double>();
        test_sort_algorithms<long double>();
        test_sort_algorithms<std::string>();
    }
    return (EXIT_SUCCESS);
}

