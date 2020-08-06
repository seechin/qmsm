#ifndef __HEAPSORTNS__
#define __HEAPSORTNS__
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
namespace HeapSortNS {
#ifndef null
  #define null NULL
#endif
    template <class HSKEY, class HSUNIT>
    class HeapSortHeap {
       public:
        int count;
        HSKEY  * e;
        HSUNIT * o;
       public:
        HeapSortHeap(){ init(0, null, null); }
        HeapSortHeap(int count, HSKEY * e, HSUNIT * o){ init(count, e, o); }
        #define Init init
        void init(int count, HSKEY * e, HSUNIT * o){ this->count = count; this->e = e; this->o = o; }
      //堆排序
       private:
        void heap_adjust(int idx, int count){  //自上向下调整一个堆
            for (int s = idx; s < count / 2; ){
                int type = 0; int k = 0;
                if (e[s] < e[2 * s + 1]) type |= 1;
                if (s < count / 2 - 1){
                    if (e[s] < e[2 * s + 2]) type |= 2;
                    if (e[2 * s + 1] < e[2 * s + 2]) type |= 4;
                }
                switch (type){
                    case 0: k = -1; break;
                    case 1: k = 2 * s + 1; break;
                    case 2: k = 2 * s + 1; break;
                    case 3: k = 2 * s + 1; break;
                    case 4: k = -1; break;
                    case 5: k = 2 * s + 2; break;
                    case 6: k = 2 * s + 2; break;
                    case 7: k = 2 * s + 2; break;
                }
                if (k >= 0){
                    HSKEY  tk = e[s]; e[s] = e[k]; e[k] = tk;
                    HSUNIT tu = o[s]; o[s] = o[k]; o[k] = tu;
                    s = k;
                }
                else return;
            }
        }
        void heap_make(int count){     //自上向下生成一个堆
            for (int i = count / 2 - 1; i >= 0; i--) heap_adjust(i, count);
        }
       public:
        #define HeapSort heap_sort
        void heap_sort(){
            heap_make(count);
            for (int i = count - 1; i>=0; ){
              //1. 把最大的记录移到后面去
                if (e[0] < e[i]) { i--; continue; }
                HSKEY  tk = e[0]; e[0] = e[i]; e[i] = tk;
                HSUNIT tu = o[0]; o[0] = o[i]; o[i] = tu;
                i--;
              //2. 对当前的堆进行调整
                heap_adjust(0,i+1);
            }
        }
    };
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#endif
