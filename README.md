# CPPSortAlgorithms
CPPSortAlgorithms is a generic C++17 header-only sorting library. In this library you can find a collection of sorting algorithms and utility, that can help you to compare speed of each one.

Sorting algorithms:

1. Bubble sort
1. Cocktail shaker sort
1. Comb sort
1. Odd–even sort
1. Insertion sort
1. Selection Sort
1. Shell sort and some variations
1. Tree sort
1. Gnome sort
1. Heap sort
1. Quick sort
1. Merge Sort
1. Bucket sort
1. Least significant digit (LSD) Radix sort
1. Most significant digit (LSD) Radix sort
1. Bitonic Sort
1. Timsort

## Build

On linux, go to project folder and type:
```
mkdir build
cd build
cmake ..
make
```

In Visual Studio open the project folder and compile the code.

After this, you can type:
```
./SortAlgorithms
```
to test each case of sorting algorithms.

Or you can type:
```
./SortAlgorithms -a comb,shell -s 10000 -t s8,u64 -sh sorted,reversed,random
```
to test specific case.

## Flags

-s -- array size;\
-a -- testing algorithms. You can type comma-separated `bubble,oddeven,shaker,comb,insertion,selection,shell,shellhib,shellpratt,shellsedgwick,shell1,shell2,shell3,tree,gnome,heap,quick,merge,bucket,lsd,msd,bitonic,tim`;\
-t -- element types in array. s8 -- signed 8 bit. u64 -- unsigned 64 bit. string -- std::string. f -- float. d -- double. ld -- long double. You can type comma-separated `s8,u8,s16,u16,s32,u32,s64,u64,f,d,ld,string`;\
-sh -- shuffle type. You can type comma separated `sorted,reversed,random,flswap,s4,s8,s16`. sorted -- sorted array. reversed -- reverse sorted array. random -- shuffled all array range. flswap -- sorted array where first and last elements are swapped (for example: 9, 2, 3, 4, 5, 6, 7, 8, 1). s4,s8,s16 -- almost sorted array where elements are shuffled in limited range (example for s4: '4, 2, 1, 3', '6, 5, 9, 8');\

Type `./SortAlgorithms -h` for more details.

## API

Bubble sort:
```
template <typename T>
sort::bubble(T start, T end);
```

Cocktail shaker sort:
```
template <typename T>
sort::shaker(T start, T end);
```

Comb sort:
```
template <typename T>
sort::comb(T start, T end);
```

Odd–even sort:
```
template <typename T>
sort::oddeven(T start, T end);
```

Insertion sort:
```
template <typename T>
sort::insertion(T start, T end);
```

Selection Sort:
```
template <typename T>
sort::selection(T start, T end);
```

Shell sort and some variations:
```
template <typename T>
sort::shell(T start, T end);
template <typename T>
sort::shellhib(T start, T end);
template <typename T>
sort::shellsedgwick(T start, T end);
template <typename T>
sort::shellpratt(T start, T end);
template <typename T>
sort::shell1(T start, T end);
template <typename T>
sort::shell2(T start, T end);
template <typename T>
sort::shell3(T start, T end);
```

Tree sort:
```
template <typename T>
sort::tree(T start, T end);
```

Gnome sort:
```
template <typename T>
sort::gnome(T start, T end);
```

Heap sort:
```
template <typename T>
sort::heap(T start, T end);
```

Quick sort:
```
template <typename T>
sort::quick(T start, T end);
```

Merge Sort:
```
template <typename T>
sort::merge(T start, T end);
```

Bucket sort:
```
template <typename T>
sort::bucket(T start, T end);
```

Least significant digit (LSD) Radix sort:
```
template <typename T>
sort::lsd(T start, T end);
```

Most significant digit (LSD) Radix sort:
```
template <typename T>
sort::msd(T start, T end);
```

Bitonic Sort:
```
template <typename T>
sort::bitonic(T start, T end);
```

Timsort:
```
template <typename T>
sort::tim(T start, T end);
```

In most algorithms you can specify a compare function in the third parameter:
```
sort::bubble<int>(it_start, it_end, [](const auto& a, const auto& b) -> bool { return (a < b); });
```

