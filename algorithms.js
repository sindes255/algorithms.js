var Algorithm = (function (){
    if(!Algorithm)return {
        REVISION: '1',
        verison: '1.0',

        /* ========================== Search algorithms ========================== */

        /* ============= Linear search ============= */

        linearSearch : function (t, A)      // t - target element,
        {                                           // A - input array.
            var n = A.length, i = 0;

            A[n] = t;

            while (A[i] !== t) i++;

            if (i < n) return i;          // In output - index of searching element.
            else return -1;               // if not, return -1.
        },

        /* ============= Binary search ============= */

        binarySearch : function (t, A)       // t - target element,
        {                                            // A - input array.
            var i = 0, j = A.length, k;

            while (i < j) {
                k = Math.floor((i + j) / 2);
                if (t <= A[k]) j = k;
                else i = k + 1;
            }

            if (A[i] === t) return i;    // In output - index of searching element.
            else return -1;              // if not, return -1.
        },

        /* ============= Substring search ============= */

        substringSearch : function (sub, str)    // sub - target substring
        {                                                 // str - input string
            var i, j, n = sub.length,
                N = str.length - n + 1;

            for (i = 0; i < N; i++) {
                j = 0;
                while (j < n && sub.charAt(j) === str.charAt(i + j)) j++;
                if (j === n) return i;
            }                                // In output - index of first substing symbol.
                                             // if not, return -1.
            return -1;
        },

        /* ========================== Sorting algorithms ========================== */

        /* ============= Bubble sort ============= */

        bubbleSort : function (A)       // A - input array, that
        {                                        // must be sorted by increase.
            var n = A.length;
            for (var i = 0; i < n - 1; i++) {
                for (var j = 0; j < n - 1 - i; j++) {
                    if (A[j + 1] < A[j]) {
                        var t = A[j + 1];
                        A[j + 1] = A[j];
                        A[j] = t;
                    }
                }
            }
            return A;    // in output - sorted by increase A array.
        },

        /* ============= Selection sort ============= */

        selectionSort : function (A)       // A - input array, that
        {                                           // must be sorted by increase.
            var n = A.length;
            for (var i = 0; i < n - 1; i++) {
                var min = i;
                for (var j = i + 1; j < n; j++) {
                    if (A[j] < A[min]) min = j;
                }
                var t = A[min];
                A[min] = A[i];
                A[i] = t;
            }
            return A;    // in output - sorted by increase A array.
        },

        /* ============= Insert sort ============= */

        insertionSort : function (A)       // A - input array, that
        {                                           // must be sorted by increase.
            var n = A.length;
            for (var i = 0; i < n; i++) {
                var v = A[i], j = i - 1;
                while (j >= 0 && A[j] > v) {
                    A[j + 1] = A[j];
                    j--;
                }
                A[j + 1] = v;
            }
            return A;    // in output - sorted by increase A array.
        },

        /* ============= Shell sort ============= */

        shellSort : function (A) {
            var n = A.length, i = Math.floor(n / 2);
            while (i > 0) {
                for (var j = 0; j < n; j++) {
                    var k = j, t = A[j];
                    while (k >= i && A[k - i] > t) {
                        A[k] = A[k - i];
                        k -= i;
                    }
                    A[k] = t;
                }
                i = (i == 2) ? 1 : Math.floor(i * 5 / 11);
            }
            return A;
        },

        /* ============= Simple Counting sort ============= */

        simpleCountingSort : function (A) {
            var n = A.length, Count = [], B = [];
            for (var i = 0; i < n; i++) Count[i] = 0;
            for (var i = 0; i < n - 1; i++) {
                for (var j = i + 1; j < n; j++) {
                    if (A[i] < A[j]) Count[j]++;
                    else Count[i]++;
                }
            }
            for (var i = 0; i < n; i++) B[Count[i]] = A[i];
            return B;
        },

        /* ============= Gap sort ============= */


        combSort : function (A)           // A - input array, that
        {                                          // must be sorted by increase.
            var n = A.length, gap = n;

            function newGap(gap)        // Supported function
            {
                gap /= 1.3;
                if (gap == 9 || gap == 10) gap = 11;
                if (gap < 1) return 1;
                return gap;
            }

            do {
                swapped = false;
                gap = newGap(gap);
                for (var i = 0; i < n - gap; ++i) {
                    if (A[i] > A[i + gap]) {
                        swapped = true;
                        var t = A[i + gap];
                        A[i + gap] = A[i];
                        A[i] = t;
                    }
                }
            } while (gap > 1 || swapped);
            return A;
        },

        /* ============= Merge sort ============= */


        mergeSort : function (A)      // A - input array, that
        {                                      // must be sorted by increase.
            function Merge(a, low, mid, high)    // Supported function
            {
                var b = new Array(high + 1 - low), h, i, j = mid + 1, k, h = low, i = 0;
                while (h <= mid && j <= high) {
                    if (a[h] <= a[j]) {
                        b[i] = a[h];
                        h++;
                    }
                    else {
                        b[i] = a[j];
                        j++;
                    }
                    i++;
                }
                if (h > mid) {
                    for (k = j; k <= high; k++) {
                        b[i] = a[k];
                        i++;
                    }
                }
                else {
                    for (k = h; k <= mid; k++) {
                        b[i] = a[k];
                        i++;
                    }
                }
                for (k = 0; k <= high - low; k++) a[k + low] = b[k];
                return a;
            }

            function merge_sort(a, low, high) {
                if (low < high) {
                    var mid = Math.floor((low + high) / 2);
                    merge_sort(a, low, mid);
                    merge_sort(a, mid + 1, high);
                    Merge(a, low, mid, high);
                }
            }

            var n = A.length;
            merge_sort(A, 0, n - 1);
            return A;
        },

        /* ============= Heap sort ============= */

        heapSort : function (A) {
            if (A.length == 0) return [];
            var n = A.length, i = Math.floor(n / 2), j, k, t;
            while (true) {
                if (i > 0) t = A[--i];
                else {
                    n--;
                    if (n == 0) return A;
                    t = A[n];
                    A[n] = A[0];
                }
                j = i;
                k = j * 2 + 1;
                while (k < n) {
                    if (k + 1 < n && A[k + 1] > A[k]) k++;
                    if (A[k] > t) {
                        A[j] = A[k];
                        j = k;
                        k = j * 2 + 1;
                    }
                    else break;
                }
                A[j] = t;
            }
        },

        /* ============= Quick sort ============= */

        quickSort : function (A) {
            if (A.length == 0) return [];
            var a = [], b = [], p = A[0];
            for (var i = 1; i < A.length; i++) {
                if (A[i] < p) a[a.length] = A[i];
                else b[b.length] = A[i];
            }
            return Algorithm.quickSort(a).concat(p, Algorithm.quickSort(b));
        },

        /* ============= Cocktail sort ============= */

        cocktailSort : function (A)    //second name is ShakerSort.
        {
            var i = 0, j = A.length - 1, s = true, t;
            while (i < j && s) {
                s = false;
                for (var k = i; k < j; k++) {
                    if (A[k] > A[k + 1]) {
                        t = A[k];
                        A[k] = A[k + 1];
                        A[k + 1] = t;
                        s = true;
                    }
                }
                j--;
                if (s) {
                    s = false;
                    for (var k = j; k > i; k--) {
                        if (A[k] < A[k - 1]) {
                            t = A[k];
                            A[k] = A[k - 1];
                            A[k - 1] = t;
                            s = true;
                        }
                    }
                }
                i++;
            }
            return A;
        },

        /* ============= Gnome sort ============= */

        gnomeSort : function (A) {
            var n = A.length, i = 1, j = 2;
            while (i < n) {
                if (A[i - 1] < A[i]) {
                    i = j;
                    j++;
                }
                else {
                    var t = A[i - 1];
                    A[i - 1] = A[i];
                    A[i] = t;
                    i--;
                    if (i == 0) {
                        i = j;
                        j++;
                    }
                }
            }
            return A;
        },

        /* ============= Natural sort ============= */

        naturalSort : function (string_array)  // A - input array with string items
        {
            var splitters = string_array.map(makeSplitter),
                sorted = splitters.sort(compareSplitters);
            return sorted.map(function (splitter) {
                return splitter.item
            });
            function makeSplitter(item) {
                return new Splitter(item)
            }

            function Splitter(item) {
                var index = 0, from = 0, parts = [], completed = false;
                this.item = item;
                var key = item;
                this.key = key;
                this.count  =  function () {
                    return parts.length;
                };
                this.part  = function (i) {
                    while (parts.length <= i && !completed) next();
                    return i < parts.length ? parts[i] : null;
                };
                function next() {
                    if (index < key.length) {
                        while (++index) {
                            var currentIsDigit = isDigit(key.charAt(index - 1)),
                                nextChar = key.charAt(index),
                                currentIsLast = index === key.length,
                                isBorder = currentIsLast || xor(currentIsDigit, isDigit(nextChar));
                            if (isBorder) {
                                var partStr = key.slice(from, index);
                                parts.push(new Part(partStr, currentIsDigit));
                                from = index;
                                break;
                            }
                        }
                    }
                    else completed = true;
                }

                function Part(text, isNumber) {
                    this.isNumber = isNumber;
                    this.value = isNumber ? Number(text) : text;
                }
            }

            function compareSplitters(sp1, sp2) {
                var i = 0;
                do {
                    var first = sp1.part(i), second = sp2.part(i);
                    if (null !== first && null !== second) {
                        if (xor(first.isNumber, second.isNumber)) {
                            return first.isNumber ? -1 : 1;
                        }
                        else {
                            var comp = compare(first.value, second.value);
                            if (comp != 0) return comp;
                        }
                    } else return compare(sp1.count(), sp2.count());
                } while (++i);
                function compare(a, b) {
                    return a < b ? -1 : a > b ? 1 : 0;
                }
            }

            function xor(a, b) {
                return a ? !b : b;
            }

            function isDigit(chr) {
                var code = charCode(chr);
                return code >= charCode("0") && code <= charCode("9");
                function charCode(c) {
                    return c.charCodeAt(0);
                }
            }
        },

        /* ========================== Searching of unique elements of an array algorithms ========================== */

        /* ============= Checking the uniqueness of elements ============= */

        testUnique : function (A) {
            var n = A.length;
            for (var i = 0; i < n - 1; i++) {
                for (var j = i + 1; j < n; j++) {
                    if (A[i] === A[j]) return false;
                }
            }
            return true;
        },

        /* ============= Removing duplicate items ============= */

        unique : function (A) {
            var n = A.length, k = 0, B = [];
            for (var i = 0; i < n; i++) {
                var j = 0;
                while (j < k && B[j] !== A[i]) j++;
                if (j == k) B[k++] = A[i];
            }
            return B;
        },

        /* ============= Getting unique items ============= */

        getUniqueElems : function (A)   // A - input ordered array.
        {
            var n = A.length, B = [];
            for (var i = 1, j = 0, t; i < n + 1; i++) {
                if (A[i - 1] === A[i]) t = A[i - 1];
                if (A[i - 1] !== t) B[j++] = A[i - 1];
            }
            // In output - array, only with
            return B;       // non-doubled elements of input array.
        },

        /* ============= Getting not unique items ============= */

        getDublsSortArr : function (A)   // A - input ordered array.
        {
            var n = A.length, B = [];

            for (var i = 1, j = 0; i < n; i++) {
                if (A[i - 1] === A[i]) B[j++] = A[i - 1];
            }

            return Algorithm.unique(B);    //   UniqueSort function see above
        },

        /* ========================== Algorithms of merge, intersection, difference and combining arrays ========================== */

        /* ============= Merging arrays ============= */

        /* ====== For two arrays ====== */

        mergeArray : function (A, B)                     // A и B - input ordered arrays.
        {
            var N = A.length, M = B.length, C = [];

            for (var i = 0, j = 0, k = 0; k < N + M; k++) {
                if (i == N) {
                    C[k] = B[j++];
                    continue;
                }
                if (j == M) {
                    C[k] = A[i++];
                    continue;
                }
                C[k] = (A[i] < B[j]) ? A[i++] : B[j++];
            }
            // In output - ordered array C,
            return C;                                 // consisting of the elements A и B.
        },

        /* ====== For a set of arrays ====== */

        multiMergeArray : function (k, A)   // Calling only with k=0. А - two-dimensional array,
        {                                           //  his elements - A[ i ] is a ordered arrays there must be merge.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.mergeArray(A[n - 2], A[n - 1]);   // merge function see above
            else
                return Algorithm.mergeArray(A[k], Algorithm.multiMergeArray(k + 1, A));  // In output - ordered one-dimensional array,
            // includes elements of input two-dimensional array.
        },


        /* ============= Intersection of arrays ============= */

        /* ====== For two arrays ====== */

        intersecArrays : function (A, B) {
            var m = A.length, n = B.length, c = 0, C = [];
            for (var i = 0; i < m; i++) {
                var j = 0, k = 0;
                while (B[j] !== A[i] && j < n) j++;
                while (C[k] !== A[i] && k < c) k++;
                if (j != n && k == c) C[c++] = A[i];
            }
            return C;
        },

        /* ====== For a set of arrays ====== */

        multiIntersecArrays : function (k, A) // Calling only with k=0. А - two-dimensional array,
        {                                             //  his elements - A[ i ] is a ordered arrays there must be crossed.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.intersecArrays(A[n - 2], A[n - 1]);   // IntersecArrays function see above
            else
                return Algorithm.intersecArrays(A[k], Algorithm.multiIntersecArrays(k + 1, A));
        },

        /* ====== For two ordered arrays ====== */

        intersecSortArr : function (A, B)         // A и B - input ordered array.
        {
            var M = A.length, N = B.length, C = [],
                m = 1, n = 1, k = 0, a = 0, b = 0;

            for (var i = 1, t = A[0]; i < M; i++)  // Shift to the beginning of the unique elements
            {
                if (A[i] !== t)                    // First m elements be unique
                {
                    A[m++] = A[i];
                    t = A[i];
                }
            }

            for (var i = 1, t = B[0]; i < N; i++)  // Similar to the previous
            {
                if (B[i] !== t) {
                    B[n++] = B[i];
                    t = B[i];
                }
            }

            while (a < m && b < n)                // Push to C not unique elements
            {
                if (A[a] < B[b]) ++a;
                else if (A[a] > B[b]) ++b;
                else С[k++] = A[a++];
            }

            return C;   // In output - C array includes A and B
                        // not unique elements, without repeats.
        },

        /* ====== For a set of ordered arrays ====== */

        multiIntersecSortArr : function (k, A)  // Calling only with k=0. А - array of two-dimensional array,
        {                                             //  his elements - A[ i ] is a ordered arrays there must be crossed.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.intersecSortArr(A[n - 2], A[n - 1]);   // intersecSortArr() function see above
            else
                return Algorithm.intersecSortArr(A[k], Algorithm.multiIntersecSortArr(k + 1, A));
        },

        /* ============= Difference of arrays ============= */

        /* ====== For two arrays ====== */

        diffArrays : function (A, B) {
            var M = A.length, N = B.length, c = 0, C = [];
            for (var i = 0; i < M; i++) {
                var j = 0, k = 0;
                while (B[j] !== A[i] && j < N) j++;
                while (C[k] !== A[i] && k < c) k++;
                if (j == N && k == c) C[c++] = A[i];
            }
            return C;
        },

        /* ====== For a set of arrays ====== */

        multiDiffArrays : function (k, A)   // Calling only with k=0. А - two-dimensional array,
        {                                           //  his elements - A[ i ] is a ordered arrays. Need to find there difference.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.diffArrays(A[n - 2], A[n - 1]);      // DiffArrays  function see above
            else
                return Algorithm.diffArrays(A[k], Algorithm.multiDiffArrays(k + 1, A));
        },                 //  In output - array A[0] with elements(without doubles),
        //  that not complicated with elements of arrays A[1],...,A[n-1].

        /* ====== For two ordered arrays ====== */

        diffSortArr : function (A, B)     // A и B - input ordered arrays.
        {
            var C = Algorithm.intersecSortArr(A, B),    // C - array of not unique
                M = A.length,                //   elements from array A and B
                N = C.length;                // IntersecSortArr()  function see above

            for (var i = 0, k = 0, a = 0, c = 0; i < M + N; i++) {
                if (A[a] === C[c]) {
                    ++a;
                    ++c;
                }
                else {
                    A[k] = A[i];
                    k++;
                    a++;
                }
            }                                        // Mathematically: removal of a plurality of subsets
            A.length = k;


            return A;     //  unique elements form A array that not includes B array
        },

        /* ====== For a set of ordered arrays ====== */

        multiDiffSortArr : function (k, A)   // Calling only with k=0. А - array of two-dimensional array,
        {                                           //  his elements - A[ i ] is a ordered arrays. Need to find there difference.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.diffSortArr(A[n - 2], A[n - 1]);      // DiffSortArr()  function see above
            else
                return Algorithm.diffSortArr(A[k], Algorithm.multiDiffSortArr(k + 1, A));
        },                 //  In output - array A[0] with elements(without doubles),
        //  that not complicated with elements of arrays A[1],...,A[n-1].

        /* ============= Symmetric difference of arrays ============= */

        /* ====== For two arrays ====== */

        symmDiffSortArr : function (A, B)     // A и B - input ordered arrays.
        {
            var N = A.length, M = B.length;

            for (var i = 1, j = 1, k = A[0]; i < N; i++)       // Delete repeating elements
            {
                if (A[i] !== k)                       // in sorted array А и B.
                {
                    A[j++] = A[i];
                    k = A[i];
                }
            }
            A.length = j;

            for (var i = 1, j = 1, k = B[0]; i < M; i++) {
                if (B[i] !== k) {
                    B[j++] = B[i];
                    k = B[i];
                }
            }
            B.length = j;

            var N = A.length, M = B.length, C = [];

            for (var i = 0, j = 0, k = 0; k < N + M; k++)        // Merge arrays A and B to C array
            {
                if (i == N) {
                    C[k] = B[j++];
                    continue;
                }
                if (j == M) {
                    C[k] = A[i++];
                    continue;
                }
                C[k] = (A[i] < B[j]) ? A[i++] : B[j++];
            }

            var N = C.length;

            for (var i = 1, j = 0, t; i < N + 1; i++)          // Shift to the beginning
            {
                if (C[i - 1] === C[i]) t = C[i - 1];        //  not- doubles elements,
                if (C[i - 1] !== t) C[j++] = C[i - 1];      //  их будет j-штук.
            }                                         // Push not-doubled j elements to C
            C.length = j;

            return C;        // In output - array C, includes elements from A and B
                             //   and not doubled in there
        },

        /* ====== For a set of ordered arrays ====== */

        multiSymmDiffSortArr : function (k, A)   // Calling only with k=0. А - array of two-dimensional array,
        {                                                //  his elements - A[ i ] is a ordered arrays. Need to find there symmetric difference.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.symmDiffSortArr(A[n - 2], A[n - 1]);     // SymmDiffSortArr()  function see above
            else
                return Algorithm.symmDiffSortArr(A[k], Algorithm.multiSymmDiffSortArr(k + 1, A));
        },

        /* ============= Combining arrays ============= */

        /* ====== For two ordered arrays ====== */

        unionSortArr : function (A, B)             // A и B - input ordered arrays.
        {
            var N = A.length, M = B.length, C = [];

            for (var i = 0, j = 0, k = 0; k < N + M; k++)        // Merge arrays A and B to C array
            {
                if (i == N) {
                    C[k] = B[j++];
                    continue;
                }
                if (j == M) {
                    C[k] = A[i++];
                    continue;
                }
                C[k] = (A[i] < B[j]) ? A[i++] : B[j++];
            }

            for (var i = 1, j = C[0], k = 1; i < N + M; i++)     // Delete repeating elements in C array
            {
                if (C[i] !== j) {
                    C[k++] = C[i];
                    j = C[i];
                }
            }
            C.length = k;
            // In output - array C, includes elements from A and B, without repeats
            return C;
        },

        /* ====== For a set of ordered arrays ====== */

        multiUnionSortArr : function (k, A)  // Calling only with k=0. А - array of two-dimensional array,
        {                                            //  his elements - A[ i ] is a ordered arrays. Need to find there symmetric difference.
            var n = A.length;
            if (k == n - 2)
                return Algorithm.unionSortArr(A[n - 2], A[n - 1]);   // UnionSortArr  function see above
            else
                return Algorithm.unionSortArr(A[k], Algorithm.multiUnionSortArr(k + 1, A));
            // In output - one-dimesional array, includes elements A[i],
        },                                      // from input matrix A without repeats.

        /* ========================== GCD & LCM(Greatest common divisor and Least common multiple) ========================== */

        /* ============= GCD ============= */

        GCD : function (A) {
            var n = A.length, x = Math.abs(A[0]);
            for (var i = 1; i < n; i++) {
                var y = Math.abs(A[i]);
                while (x && y) {
                    x > y ? x %= y : y %= x;
                }
                x += y;
            }
            return x;
        },

        /* ============= LCM ============= */

        LCM : function (A) {
            var n = A.length, a = Math.abs(A[0]);
            for (var i = 1; i < n; i++) {
                var b = Math.abs(A[i]), c = a;
                while (a && b) {
                    a > b ? a %= b : b %= a;
                }
                a = Math.abs(c * A[i]) / (a + b);
            }
            return a;
        },

        NOD : function (A) {
            return Algorithm.GCD(A)
        },

        NOK : function (A) {
            return Algorithm.LCM(A)
        },

        /* ========================== Operations on two-dimensional matrix) ========================== */

        /* ============= Matrix transposition ============= */

        transMatrix : function (A)       //two-dimensional input matrix
        {
            var m = A.length, n = A[0].length, AT = [];
            for (var i = 0; i < n; i++) {
                AT[i] = [];
                for (var j = 0; j < m; j++) AT[i][j] = A[j][i];
            }
            return AT;
        },

        /* ============= Addition of matrices ============= */

        sumMatrix : function (A, B)       //two two-dimensional input matrices
        {
            var m = A.length, n = A[0].length, C = [];
            for (var i = 0; i < m; i++) {
                C[i] = [];
                for (var j = 0; j < n; j++) C[i][j] = A[i][j] + B[i][j];
            }
            return C;
        },

        /* ============= Multiplication between matrices and numbers ============= */

        multMatrixNumber : function (a, A)  // a - number, A - two-dimensional input matrix
        {
            var m = A.length, n = A[0].length, B = [];
            for (var i = 0; i < m; i++) {
                B[i] = [];
                for (var j = 0; j < n; j++) B[i][j] = a * A[i][j];
            }
            return B;
        },

        /* ============= Multiplication of matrices ============= */

        multiplyMatrix : function (A, B) {
            var rowsA = A.length, colsA = A[0].length,
                rowsB = B.length, colsB = B[0].length,
                C = [];
            if (colsA != rowsB) return false;
            for (var i = 0; i < rowsA; i++) C[i] = [];
            for (var k = 0; k < colsB; k++) {
                for (var i = 0; i < rowsA; i++) {
                    var t = 0;
                    for (var j = 0; j < rowsB; j++) t += A[i][j] * B[j][k];
                    C[i][k] = t;
                }
            }
            return C;
        },

        /* ============= Powering of matrices ============= */

        matrixPow : function (n, A) {
            if (n == 1) return A;     // MultiplyMatrix function see above
            else return Algorithm.multiplyMatrix(A, Algorithm.matrixPow(n - 1, A));
        },

        /* ============= Matrix determinant ============= */

        determinantMatrix : function (A)   // Using Bareys's algorithm, difficult O(n^3)
        {
            var N = A.length, B = [], denom = 1, exchanges = 0;
            for (var i = 0; i < N; ++i) {
                B[i] = [];
                for (var j = 0; j < N; ++j) B[i][j] = A[i][j];
            }
            for (var i = 0; i < N - 1; ++i) {
                var maxN = i, maxValue = Math.abs(B[i][i]);
                for (var j = i + 1; j < N; ++j) {
                    var value = Math.abs(B[j][i]);
                    if (value > maxValue) {
                        maxN = j;
                        maxValue = value;
                    }
                }
                if (maxN > i) {
                    var temp = B[i];
                    B[i] = B[maxN];
                    B[maxN] = temp;
                    ++exchanges;
                }
                else {
                    if (maxValue == 0) return maxValue;
                }
                var value1 = B[i][i];
                for (var j = i + 1; j < N; ++j) {
                    var value2 = B[j][i];
                    B[j][i] = 0;
                    for (var k = i + 1; k < N; ++k) B[j][k] = (B[j][k] * value1 - B[i][k] * value2) / denom;
                }
                denom = value1;
            }
            if (exchanges % 2) return -B[N - 1][N - 1];
            else return B[N - 1][N - 1];
        },

        /* ============= Matrix rank ============= */

        matrixRank : function (A) {
            var m = A.length, n = A[0].length, k = (m < n ? m : n), r = 1, rank = 0;
            while (r <= k) {
                var B = [];
                for (var i = 0; i < r; i++) B[i] = [];
                for (var a = 0; a < m - r + 1; a++) {
                    for (var b = 0; b < n - r + 1; b++) {
                        for (var c = 0; c < r; c++) {
                            for (var d = 0; d < r; d++) B[c][d] = A[a + c][b + d];
                        }
                        if (Algorithm.determinantMatrix(B) != 0) rank = r;
                    }       // determinantMatrix  function see above
                }
                r++;
            }
            return rank;
        },

        /* ============= Adjugate matrix ============= */

        adjugateMatrix : function (A)   // A - two-dimensional input matrix
        {
            var N = A.length, adjA = [];
            for (var i = 0; i < N; i++) {
                adjA[i] = [];
                for (var j = 0; j < N; j++) {
                    var B = [], sign = ((i + j) % 2 == 0) ? 1 : -1;
                    for (var m = 0; m < j; m++) {
                        B[m] = [];
                        for (var n = 0; n < i; n++)   B[m][n] = A[m][n];
                        for (var n = i + 1; n < N; n++) B[m][n - 1] = A[m][n];
                    }
                    for (var m = j + 1; m < N; m++) {
                        B[m - 1] = [];
                        for (var n = 0; n < i; n++)   B[m - 1][n] = A[m][n];
                        for (var n = i + 1; n < N; n++) B[m - 1][n - 1] = A[m][n];
                    }
                    adjA[i][j] = sign * Algorithm.determinantMatrix(B);   // Determinant function see above
                }
            }
            return adjA;
        },

        /* ============= Inverse matrix ============= */

        inverseMatrix : function (A)   // A - two-dimensional input matrix
        {
            var det = Algorithm.determinantMatrix(A);                // Determinant function see above
            if (det == 0) return false;
            var N = A.length, A = Algorithm.adjugateMatrix(A); //AdjugateMatrix function see above
            for (var i = 0; i < N; i++) {
                for (var j = 0; j < N; j++) A[i][j] /= det;
            }
            return A;
        },


        /* ========================== Randomizer algorithms ========================== */

        /* ============= Random integer between min and max ============= */

        getRandomArbitary : function (min, max) {
            return Math.random() * (max - min) + min;
        },


        /* ============= Random number between min and max ============= */


// uses Math.round() gives maldistribution!
        getRandomInt : function (min, max) {
            return Math.floor(Math.random() * (max - min + 1)) + min;
        },

        /* ========================== Math algorithms ========================== */

        /* ============= Get number from Fibonacci by his position ============= */

        getFibonahhiNumber : function (n) {// n - number of position
            var result;
            //using Bine formula
            result = (Math.pow(((1 + Math.sqrt(5)) / 2), n) - Math.pow(((1 - Math.sqrt(5)) / 2), n)) / Math.sqrt(5)

            return result
        },


        /* ============= Get position from Fibonacci by number  ============= */

        getFibonahhiPosition : function (n) {// n - number
            var result, a, b, c, mod;
            if (n < 3 && n > -3) {
                result = n;
            } else {
                a = 1;
                b = 1;
                c = a + b;
                mod = 1;
                if (n < 0)mod = -1;
                n = Math.abs(n);
                for (var i = 2; c < n; i++) {
                    c = a + b;
                    a = b;
                    b = c;
                }
                if (c > n) {
                    result = false
                } else {
                    result = i * mod;
                }
            }
            return result;//return false if not Fibonahhi number
        }
    }
}());







