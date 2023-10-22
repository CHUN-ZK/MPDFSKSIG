#pragma once

#define NSYM 0
#define SYM 1

/*��Xdetail*/inline int verbose = TRUE; /* default -v: ��X�ثe���窺detail */

/* ����*/
enum keys {
    KEY_SYMMETRY, KEY_KSIG, KEY_TEST_SYM, KEY_TEST_KSIG, KEY_TEST_LARGECIRCUIT, KEY_CHUNK, KEY_THREAD, KEY_ORDERING, KEY_SYMRATE
};
struct {
    string name;
    enum keys key;
} option_table[] = {
    {"symmetry", KEY_SYMMETRY},        /* ��ٰt���ˬd */
    {"ksig", KEY_KSIG},                /* K-�S�x�ȭp�� */
    {"testsym", KEY_TEST_SYM},         /* ��MCNC�BISCAS���չq���i�����ˬd���� */
    {"testLC", KEY_TEST_LARGECIRCUIT}, /* ��ͦ����j���q���i�����ˬd���� */
    {"testksig", KEY_TEST_KSIG},       /* ��MCNC���չq���i��K-�S�x�ȹ��� */
    {"testchunk", KEY_CHUNK},           /* �ϥΤ��Pchunksize�i��DCP����ˬd���� */
    {"testthread", KEY_THREAD},          /* �ϥΤ��P������ƶq�i��DCP����ˬd����*/
    {"testorder", KEY_ORDERING},        /* �ϥΤ��P��ƶ���J�u�@��C���覡�i��DCP����ˬd����*/
    {"testrate", KEY_SYMRATE},         /* �ϥγ�Ӱ�����̧ǰ����ƶ�������ٰt�諸��ҹ���*/
};

/* ��@��k */
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