# MPDFSKSIG

布林函數對稱性偵測與K-特徵值計算的多執行緒編程(Multithreaded Programming for Detectin Funtional Symmetries and Computing K-Signatures)

# 使用說明

-D 使用子命令

  sym 對輸入或生成電路進行對稱檢查
  
  ksig 對輸入電路進行K-特徵值計算，生成電路不可使用
  
  testsym 使用Naive、Intel TBB、OpenMP與DCP(Ours)對電路進行對稱檢查並輸出執行時間
  
  testthread 在使用DCP進行對稱檢查時，檢查1、2、4、8、16與24個最大使用執行緒對執行時間的影響
  
  testchunk 在使用DCP進行對稱檢查時，檢查多種chunksize對執行時間的影響
  
  testorder 在使用DCP進行對稱檢查時，檢查以橫向、斜向與權重的方式將資料塊放入工作佇列種對執行時間的影響
  
  testksig 使用Naive、Intel TBB、OpenMP與DCP(Ours)對電路進行K-特徵值計算並輸出執行時間

-R 輸入PLA格式更改為FR_type，在輸入生成的電路時所使用

-M 選擇實作的方法，在使用子命令sym與ksig時使用

  naive 使用Naive的實作方法
  tbb 使用Intel TBB的實作方法
  openMP 使用openMP的實作方法
  dcp 使用DCP的實作方法

-t 設定最大的執行緒數量，子命令testthread不會使用此參數

  ex: -t 20 -> 最大執行緒數量為20

-c 設定chunksize的大小，子命令testchunk不會使用此參數
  ex: -c 500 500 -> chunksize為500 * 500

-o 選擇使用的資料塊放入工作佇列方式

  h 使用橫向(HORIZONTAL)的方式
  
  d 使用斜向(DIAGONALLY)的方式
  
  w 使用權重(WEIGHT)的方式

-K 選擇K-特徵值(KSO或KSI)

  kso 計算KSO
  
  ksi 計算KSI

-G 選擇使用生成電路進行實驗，計算K-特徵值的實驗皆不能使用

-S 生成沒有函數對稱性的電路。如果沒有輸入此參數，那將會生成有函數對稱性的電路

-i 生成電路的輸入數量

  ex: -i 1000 -> 輸入數量為1000
  
-f 生成電路的on-set積項數量

  ex: -f 5000 -> on-set積項數量為5000
  
-r 生成電路的off-set積項數量

  ex: -r 5000 -> off-set積項數量為5000
  
-h 生成電路中HD<3的積項配對數量

  ex: -h 5000 -> HD<3的積項配對數量為5000

-v 關閉輸出的detail。如果沒有輸入此參數，那輸出實驗的detail

注意: 一定要輸入pla檔，且一定要放在參數的最後一位，在生成電路的情況下，將會把生成的電路寫入到輸入的pla檔中

# 使用Library

## Espresso
Library來源為: https://github.com/Gigantua/Espresso

## Intel TBB
Library來源為: https://github.com/oneapi-src/oneTBB ，
使用的版本為oneapi-tbb-2021.9.0

## boost
需自行安裝並配置，網址為 https://www.boost.org/ ，
使用的版本為boost_1_82_0

# 聯絡
若有疑問，請寄信至f410226036@gmail.com
