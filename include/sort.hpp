/*!
 * @file       sort.hpp
 * @brief      Defines the sorting algorithms
 * @author     Alexander Alexeev &lt;ntllct@protonmail.com&gt;
 * @date       Oct 15, 2018
 * @copyright  Copyright &copy; 2018 Alexander Alexeev.
 *             This project is released under the MIT License.
 */

/*********************************************************************
 * MIT License
 *
 * Copyright (c) 2018 Alexander Alexeev [ntllct@protonmail.com]
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *********************************************************************/

#pragma once

#include <cmath>
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <stack>
#include <type_traits>
#include <vector>

namespace sort {

template <class T>
class Heap {
 public:
  size_t size() { return (nElements); }
  T top() { return (h[0]); }
  bool empty() { return (nElements == 0); }
  void push(T a) {
    h.push_back(a);
    up(nElements++);
  }
  void pop() {
    std::swap(h[--nElements], h[0]);
    h.pop_back();
    down(0);
  }
  void clear() {
    h.clear();
    nElements = 0;
  }
  void reserve(size_t size) { h.reserve(size); }
  void shrink_to_fit() { h.shrink_to_fit(); }
  T operator[](size_t a) { return h[a]; }
  Heap() = default;
  Heap(size_t size) { h.reserve(size); }
  Heap(const Heap&) = default;
  Heap(Heap&&) = default;
  Heap& operator=(const Heap&) = default;
  Heap& operator=(Heap&&) = default;
  ~Heap() = default;

 private:
  std::vector<T> h{};
  size_t nElements{0};

  void up(size_t a) {
    while (a > 0) {
      size_t p = (a - 1) / 2;
      if (h[p] > h[a])
        std::swap(h[p], h[a]);
      else
        break;
      a = p;
    }
  }
  void down(size_t a) {
    size_t next_level = 2 * a + 1;
    while (next_level < nElements) {
      size_t l = next_level;
      size_t r = next_level + 1;
      if (r == nElements) {
        if (h[l] < h[a]) std::swap(h[l], h[a]);
        break;
      } else if (h[l] <= h[r]) {
        if (h[l] < h[a]) {
          std::swap(h[l], h[a]);
          a = l;
        } else {
          break;
        }
      } else if (h[r] < h[a]) {
        std::swap(h[r], h[a]);
        a = r;
      } else {
        break;
      }
      next_level = 2 * a + 1;
    }
  }
};

auto default_compare = [](const auto& a, const auto& b) -> bool {
  return (a < b);
};
using COMPARE = decltype(default_compare);

template <typename T>
void bubble(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  bool need_sort = true;
  while (need_sort) {
    need_sort = false;
    for (decltype(length) i = 1; i < length; ++i) {
      if (cmps(start[i], start[i - 1])) {
        std::swap(start[i], start[i - 1]);
        need_sort = true;
      }
    }
    length--;
  }
}
template <typename T>
void oddeven(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  bool need_sort = true;
  while (need_sort) {
    need_sort = false;
    for (decltype(length) j = 1; j < length; j += 2) {
      if (cmps(start[j], start[j - 1])) {
        std::swap(start[j], start[j - 1]);
        need_sort = true;
      }
    }
    for (decltype(length) j = 2; j < length; j += 2) {
      if (cmps(start[j], start[j - 1])) {
        std::swap(start[j], start[j - 1]);
        need_sort = true;
      }
    }
  }
}
template <typename T>
void shaker(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  bool need_sort = true;
  decltype(length) start_id = 1;
  decltype(length) end_id = length;
  while (need_sort) {
    need_sort = false;
    for (decltype(length) i = start_id; i < end_id; ++i) {
      if (cmps(start[i], start[i - 1])) {
        std::swap(start[i], start[i - 1]);
        need_sort = true;
      }
    }
    if (!need_sort) break;
    end_id--;

    need_sort = false;
    for (decltype(length) i = end_id - 1; i >= start_id; --i) {
      if (cmps(start[i], start[i - 1])) {
        std::swap(start[i], start[i - 1]);
        need_sort = true;
      }
    }
    start_id++;
  }
}
template <typename T>
void comb(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  auto step = length - 1;
  float k = 1.24733F;  // Comb factor
  float fstep = static_cast<float>(step);
  while (step > 1) {
    T start_it = start;
    T start_step_it = start + step;
    while (start_step_it != end) {
      if (cmps(*start_step_it, *start_it)) std::swap(*start_it, *start_step_it);
      start_it++;
      start_step_it++;
    }
    fstep /= k;
    step = static_cast<decltype(length)>(std::floor(fstep));
  }
  bubble<T>(start, end, cmps);
}
// The only one step of a bubble sort
template <typename T>
void bubble_up(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  for (decltype(length) i = length - 1; i != 0; --i)
    if (cmps(start[i], start[i - 1])) std::swap(start[i], start[i - 1]);
}
template <typename T>
void insertion(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  for (decltype(length) i = 1; i < length; ++i)
    for (decltype(length) j = i; j != 0 && cmps(start[j], start[j - 1]); --j)
      std::swap(start[j], start[j - 1]);
}
template <typename T>
void shell(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  for (decltype(length) step = length / 2; step != 0; step >>= 1)
    for (decltype(length) i = step; i < length; ++i)
      for (decltype(length) j = i; j >= step && cmps(start[j], start[j - step]);
           j -= step)
        std::swap(start[j], start[j - step]);
}
template <typename T>
void shellhib(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  decltype(length) step = 1;
  while (step < length) step <<= 1;
  step >>= 1;
  if (step > 1) step--;
  for (; step != 0; step >>= 1)
    for (decltype(length) i = step; i < length; ++i)
      for (decltype(length) j = i; j >= step && cmps(start[j], start[j - step]);
           j -= step)
        std::swap(start[j], start[j - step]);
}
template <typename T>
void shellsedgwick(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  decltype(length) steps[100];
  steps[0] = 1;
  long q = 1;
  while (steps[q - 1] * 3 < length) {
    if (q % 2 == 0)
      steps[q] = 9 * (1UL << q) - 9 * (1UL << (q / 2)) + 1;
    else
      steps[q] = 8 * (1UL << q) - 6 * (1UL << ((q + 1) / 2)) + 1;
    q++;
  }
  q--;

  for (; q >= 0; --q)
    for (decltype(length) i = steps[q]; i < length; ++i)
      for (decltype(length) j = i;
           j >= steps[q] && cmps(start[j], start[j - steps[q]]); j -= steps[q])
        std::swap(start[j], start[j - steps[q]]);
}
template <typename T>
void shellpratt(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  decltype(length) steps[100];
  steps[0] = 1;
  long q = 1;
  for (decltype(length) i = 1; i < length; ++i) {
    decltype(length) curr = 1UL << i;
    if (curr > length / 2) break;
    for (decltype(length) j = 1; j < length; ++j) {
      curr *= 3;
      if (curr > length / 2) break;
      steps[q++] = curr;
    }
  }
  insertion(steps, steps + q);
  q--;

  for (; q >= 0; --q)
    for (decltype(length) i = steps[q]; i < length; ++i)
      for (decltype(length) j = i;
           j >= steps[q] && cmps(start[j], start[j - steps[q]]); j -= steps[q])
        std::swap(start[j], start[j - steps[q]]);
}
template <typename T>
void shell1(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  decltype(length) steps[100];
  steps[0] = 1;
  int q = 1;
  while (steps[q - 1] < length) {
    auto s = steps[q - 1];
    steps[q++] = s * 4 + s / 4;
  }
  q -= 2;

  for (; q >= 0; --q)
    for (decltype(length) i = steps[q]; i < length; ++i)
      for (decltype(length) j = i;
           j >= steps[q] && cmps(start[j], start[j - steps[q]]); j -= steps[q])
        std::swap(start[j], start[j - steps[q]]);
}
template <typename T>
void shell2(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  decltype(length) steps[100];
  steps[0] = 1;
  long q = 1;
  while (steps[q - 1] < length) {
    auto s = steps[q - 1];
    steps[q++] = s * 3 + s / 3;
  }
  q -= 2;

  for (; q >= 0; --q)
    for (decltype(length) i = steps[q]; i < length; ++i)
      for (decltype(length) j = i;
           j >= steps[q] && cmps(start[j], start[j - steps[q]]); j -= steps[q])
        std::swap(start[j], start[j - steps[q]]);
}
template <typename T>
void shell3(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  decltype(length) steps[100];
  steps[0] = 1;
  long q = 1;
  while (steps[q - 1] < length) {
    auto s = steps[q - 1];
    steps[q++] = s * 4 - s / 5;
  }
  q -= 2;

  for (; q >= 0; --q)
    for (decltype(length) i = steps[q]; i < length; ++i)
      for (decltype(length) j = i;
           j >= steps[q] && cmps(start[j], start[j - steps[q]]); j -= steps[q])
        std::swap(start[j], start[j - steps[q]]);
}
template <typename T>
void tree(T start, T end, bool is_reversed = false) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  std::multiset<decltype(typename std::iterator_traits<T>::value_type())> tree{};
  for (auto it = start; it != end; ++it) tree.insert(*it);
  if (!is_reversed) {
    for (auto it = tree.cbegin(); it != tree.cend(); ++it) *start++ = *it;
  } else {
    for (auto it = tree.crbegin(); it != tree.crend(); ++it) *start++ = *it;
  }
}
template <typename T>
void gnome(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  T it = start + 1;

  while (it != end) {
    if (it == start || !cmps(*it, *(it - 1))) {
      it++;
    } else {
      std::swap(*(it - 1), *it);
      it--;
    }
  }
}
template <typename T>
void selection(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  decltype(typename std::iterator_traits<T>::value_type()) val{};
  for (decltype(length) i = 0; i < length; ++i) {
    val = start[i];
    decltype(length) pos = i;
    for (decltype(length) j = i + 1; j < length; ++j)
      if (cmps(start[j], val)) {
        val = start[j];
        pos = j;
      }
    start[pos] = start[i];
    start[i] = val;
  }
}
template <typename T>
void heap(T start, T end) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  Heap<decltype(typename std::iterator_traits<T>::value_type())> h(length);
  T it = start;
  while (it != end) h.push(*it++);
  it = start;
  while (it != end) {
    *it++ = h.top();
    h.pop();
  }
}
template <typename T>
void quick(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  decltype(typename std::iterator_traits<T>::value_type()) support_val =
      *(start + length / 2);
  T left_it = start;
  T right_it = end - 1;
  while (left_it < right_it) {
    while (cmps(*left_it, support_val)) left_it++;
    while (cmps(support_val, *right_it)) right_it--;
    if (left_it <= right_it) {
      std::swap(*left_it, *right_it);
      left_it++;
      right_it--;
    }
  }
  if (start < right_it) quick(start, right_it + 1, cmps);
  if (left_it < end) quick(left_it, end, cmps);
}
template <typename T>
void merge(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  using TN = decltype(typename std::iterator_traits<T>::value_type());
  auto temp = std::make_unique<TN[]>(length);

  std::function<void(T, T, COMPARE)> fn_mergesort = nullptr;
  fn_mergesort = [&fn_mergesort, &temp](T start, T end,
                                        COMPARE cmps = default_compare) {
    if (start >= end) return;
    auto length = end - start;
    if (length < 2) return;
    T middle = start + length / 2;
    fn_mergesort(start, middle, cmps);
    fn_mergesort(middle, end, cmps);
    std::merge(start, middle, middle, end, temp.get());
    std::copy_n(temp.get(), std::distance(start, end), start);
  };
  fn_mergesort(start, end, cmps);
}
// fix float
template <typename T>
typename std::enable_if<
    std::is_arithmetic<
        decltype(typename std::iterator_traits<T>::value_type())>::value,
    void>::type
bucket(T start, T end, COMPARE cmps = default_compare) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  using TN = decltype(typename std::iterator_traits<T>::value_type());
  std::unique_ptr<TN[]> temp = std::make_unique<TN[]>(length);

  std::function<void(T, T, COMPARE)> fn_busket = nullptr;
  fn_busket = [&fn_busket, &temp](T start, T end, COMPARE cmps) {
    if (end - start <= 64) {
      insertion(start, end, cmps);
      return;
    }
    // Find min and max elements, and check sorted
    TN minz = *start;
    TN maxz = *start;
    bool is_sorted = true;
    for (T it = start + 1; it < end; ++it) {
      minz = std::min(minz, *it);
      maxz = std::max(maxz, *it);
      if (cmps(*it, *(it - 1))) is_sorted = false;
    }
    if (is_sorted) return;
    auto diff = maxz - minz + 1;
    size_t numbuckets = (end - start <= 1e7) ? 1500 : 3000;
    auto range = (diff + numbuckets - 1) / numbuckets;
    auto cnt = std::make_unique<size_t[]>(numbuckets + 1);
    std::fill_n(cnt.get(), numbuckets + 1, 0);
    size_t cur = 0;
    for (T i = start; i < end; ++i) {
      temp[cur++] = *i;
      auto ind = (*i - minz) / range;
      cnt[ind + 1]++;
    }
    size_t sz = 0;
    for (size_t i = 1; i <= numbuckets; ++i)
      if (cnt[i]) sz++;
    auto run = std::make_unique<size_t[]>(sz);
    cur = 0;
    for (size_t i = 1; i <= numbuckets; ++i)
      if (cnt[i]) run[cur++] = i - 1;
    for (size_t i = 1; i <= numbuckets; ++i) cnt[i] += cnt[i - 1];
    cur = 0;
    for (T i = start; i < end; ++i) {
      auto ind = (temp[cur] - minz) / range;
      *(start + cnt[ind]) = temp[cur];
      cur++;
      cnt[ind]++;
    }
    for (size_t i = 0; i < sz; ++i) {
      auto r = run[i];
      if (r != 0)
        fn_busket(start + cnt[r - 1], start + cnt[r], cmps);
      else
        fn_busket(start, start + cnt[r], cmps);
    }
  };
  fn_busket(start, end, cmps);
}
// Least Significant Digit
template <typename T>
typename std::enable_if<std::is_unsigned<decltype(typename std::iterator_traits<
                                                  T>::value_type())>::value,
                        void>::type
lsd(T start, T end) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  using TN = decltype(typename std::iterator_traits<T>::value_type());
  const size_t N = 8;
  const int MASK = 0xff;
  auto digit = [N, MASK](TN n, int k) -> TN { return (n >> (N * k) & MASK); };

  auto temp = std::make_unique<TN[]>(length);
  auto counters = std::make_unique<size_t[]>(MASK + 1);
  for (unsigned int i = 0; i < sizeof(TN); ++i) {
    std::fill_n(counters.get(), MASK + 1, 0);
    for (auto it = start; it < end; ++it) counters[digit(*it, i)]++;
    for (int j = 1; j < MASK + 1; ++j) counters[j] += counters[j - 1];
    for (auto j = end - 1; j >= start; --j) temp[--counters[digit(*j, i)]] = *j;
    std::copy_n(temp.get(), length, start);
  }
}
// Most Significant Digit
template <typename T>
typename std::enable_if<std::is_unsigned<decltype(typename std::iterator_traits<
                                                  T>::value_type())>::value,
                        void>::type
msd(T start, T end) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  if (length <= 32) {
    insertion(start, end);
    return;
  }

  using TN = decltype(typename std::iterator_traits<T>::value_type());
  auto temp = std::make_unique<TN[]>(length);
  const size_t N = 8;
  const size_t MASK = 0xff;
  auto digit = [N, MASK](TN n, int k) -> TN { return (n >> (N * k) & MASK); };
  std::function<void(T, T, int)> radix_msd = nullptr;
  radix_msd = [&digit, &radix_msd, MASK, &temp](T start, T end, int d) {
    if (d == -1) return;
    auto length = std::distance(start, end);
    if (length <= 32) {
      insertion(start, end);
      return;
    }
    auto cnt = std::make_unique<size_t[]>(MASK + 2);
    std::fill_n(cnt.get(), MASK + 2, 0);
    std::copy(start, end, temp.get());
    for (T it = start; it < end; ++it) cnt[digit(*it, d) + 1]++;
    size_t sz = std::count_if(cnt.get() + 1, cnt.get() + MASK + 2,
                              [](auto val) { return (val != 0); });
    auto run = std::make_unique<size_t[]>(sz);
    size_t cur = 0;
    for (size_t i = 1; i <= MASK + 1; ++i)
      if (cnt[i] != 0) run[cur++] = i - 1;
    for (size_t i = 1; i <= MASK + 1; ++i) cnt[i] += cnt[i - 1];
    cur = 0;
    for (T i = start; i < end; ++i, ++cur) {
      size_t ind = digit(temp[cur], d);
      *(start + cnt[ind]++) = temp[cur];
      if (start + cnt[ind] > end) std::cout << "!!!!!!!" << std::endl;
    }
    for (size_t i = 0, r = run[i]; i < sz; ++i) {
      r = run[i];
      auto ss = (r != 0) ? start + cnt[r - 1] : start;
      radix_msd(ss, start + cnt[r], d - 1);
    }
  };
  radix_msd(start, end, sizeof(TN) - 1);
}

template <typename T>
typename std::enable_if<std::is_unsigned<decltype(typename std::iterator_traits<
                                                  T>::value_type())>::value,
                        void>::type
bitonic(T start, T end) {
  auto length = std::distance(start, end);
  if (length < 2) return;

  using TN = decltype(typename std::iterator_traits<T>::value_type());

  std::function<void(TN*, TN*, bool)> bitseq = nullptr;
  bitseq = [&bitseq](TN* start, TN* end, bool inv) -> void {
    auto length = end - start;
    if (length < 2) return;
    TN* medium = start + length / 2;
    TN* it1 = start;
    TN* it2 = medium;
    while (it1 < medium && it2 < end) {
      if (inv ^ (*it1 > *it2)) std::swap(*it1, *it2);
      it1++;
      it2++;
    }
    bitseq(start, medium, inv);
    bitseq(medium, end, inv);
  };

  std::function<void(TN*, TN*)> makebitonic = nullptr;
  makebitonic = [&makebitonic, &bitseq](TN* start, TN* end) {
    auto length = end - start;
    if (length < 2) return;
    TN* medium = start + length / 2;
    makebitonic(start, medium);
    bitseq(start, medium, false);
    makebitonic(medium, end);
    bitseq(medium, end, true);
  };

  using TT = typename std::iterator_traits<T>::difference_type;

  TT n = 1;
  while (n < length) n <<= 1;
  auto temp = std::make_unique<TN[]>(n);
  auto it = std::copy(start, end, temp.get());
  TN inf = *std::max_element(start, end) + 1;
  std::fill(it, temp.get() + n, inf);
  makebitonic(temp.get(), temp.get() + n);
  bitseq(temp.get(), temp.get() + n, false);
  std::copy_n(temp.get(), length, start);
}
template <typename T>
void tim(T start, T end) {
  auto length = std::distance(start, end);
  if (length < 2) return;
  if (length <= 64) {
    insertion(start, end);
    return;
  }
  using TN = decltype(typename std::iterator_traits<T>::value_type());
  auto temp = std::make_unique<TN[]>(length);

  auto minrun = length;
  size_t f = 0;
  while (minrun >= 64) {
    f |= minrun & 1UL;
    minrun >>= 1;
  }
  minrun += f;

  T array_it = start;
  std::stack<std::pair<size_t, T>> stack_it;
  while (array_it < end) {
    T c1 = array_it;
    while (c1 < end - 1 && *c1 <= *(c1 + 1)) c1++;
    T c2 = array_it;
    while (c2 < end - 1 && *c2 >= *(c2 + 1)) c2++;
    if (c1 >= c2) {
      c1 = std::max(c1, array_it + minrun - 1);
      c1 = std::min(c1, end - 1);
      insertion(array_it, c1 + 1);
      stack_it.emplace(c1 - array_it + 1, array_it);
      array_it = c1 + 1;
    } else {
      c2 = std::max(c2, array_it + minrun - 1);
      c2 = std::min(c2, end - 1);
      std::reverse(array_it, c2 + 1);
      insertion(array_it, c2 + 1);
      stack_it.emplace(c2 - array_it + 1, array_it);
      array_it = c2 + 1;
    }
    while (stack_it.size() > 2) {
      auto x = stack_it.top();
      stack_it.pop();
      auto y = stack_it.top();
      stack_it.pop();
      auto z = stack_it.top();
      stack_it.pop();
      if (z.first >= x.first + y.first && y.first >= x.first) {
        stack_it.emplace(z);
        stack_it.emplace(y);
        stack_it.emplace(x);
        break;
      } else if (z.first >= x.first + y.first) {
        std::merge(y.second, x.second, x.second, x.second + x.first,
                   temp.get());
        std::copy_n(temp.get(), y.first + x.first, y.second);
        stack_it.emplace(z);
        stack_it.emplace(x.first + y.first, y.second);
      } else {
        std::merge(z.second, y.second, y.second, y.second + y.first,
                   temp.get());
        std::copy_n(temp.get(), z.first + y.first, z.second);
        stack_it.emplace(z.first + y.first, z.second);
        stack_it.emplace(x);
      }
    }
  }
  while (stack_it.size() > 1) {
    auto x = stack_it.top();
    stack_it.pop();
    auto y = stack_it.top();
    stack_it.pop();
    auto dst = (x.second < y.second) ? x.second : y.second;
    std::merge(y.second, y.second + y.first, x.second, x.second + x.first,
               temp.get());
    std::copy_n(temp.get(), y.first + x.first, dst);
    stack_it.emplace(y.first + x.first, dst);
  }
}
}  // namespace sort
