#pragma once

#define NSYM 0
#define SYM 1

/*輸出detail*/inline int verbose = TRUE; /* default -v: 輸出目前實驗的detail */

/* 實驗*/
enum keys {
    KEY_SYMMETRY, KEY_KSIG, KEY_TEST_SYM, KEY_TEST_KSIG, KEY_TEST_LARGECIRCUIT, KEY_CHUNK, KEY_THREAD, KEY_ORDERING, KEY_SYMRATE
};
struct {
    string name;
    enum keys key;
} option_table[] = {
    {"symmetry", KEY_SYMMETRY},        /* 對稱配對檢查 */
    {"ksig", KEY_KSIG},                /* K-特徵值計算 */
    {"testsym", KEY_TEST_SYM},         /* 對MCNC、ISCAS測試電路進行對稱檢查實驗 */
    {"testLC", KEY_TEST_LARGECIRCUIT}, /* 對生成的大型電路進行對稱檢查實驗 */
    {"testksig", KEY_TEST_KSIG},       /* 對MCNC測試電路進行K-特徵值實驗 */
    {"testchunk", KEY_CHUNK},           /* 使用不同chunksize進行DCP對稱檢查實驗 */
    {"testthread", KEY_THREAD},          /* 使用不同執行緒數量進行DCP對稱檢查實驗*/
    {"testorder", KEY_ORDERING},        /* 使用不同資料塊放入工作佇列的方式進行DCP對稱檢查實驗*/
    {"testrate", KEY_SYMRATE},         /* 使用單個執行緒依序執行資料塊移除對稱配對的比例實驗*/
};

/* 實作方法 */
enum implements {
    ALL, NAIVE, TBB, OPENMP, DCP, WRONG
};
struct {
    string name;
    enum implements key;
    int id;
    string full_name;
} implements_table[] = {
    {"all", ALL, -1, "ALL"},
    {"naive", NAIVE, 0, "NAIVE"},
    {"tbb", TBB, 2, "TBB"},
    {"openMP", OPENMP, 3, "OPENMP"},
    {"dcp", DCP, 1, "DCP"}
};

enum ORDER {
    HORIZONTAL, DIAGONALLY, WEIGHT
};
struct {
    string name;
    enum ORDER method;
    string full_name;
} ordering_table[] = {
    {"h", HORIZONTAL, "HORIZONTAL"},
    {"d", DIAGONALLY, "DIAGONALLY"},
    {"w", WEIGHT, "WEIGHT"}
};

struct {
    string name;
    int k_mode;
    string full_name;
} ksig_table[] = {
    {"kso", 0, "KSO"},
    {"ksi", 1, "KSI"},
    {"all", 2, "KSO and KSI"},
};