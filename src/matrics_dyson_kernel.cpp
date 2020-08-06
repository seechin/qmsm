//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Simple Matrix Algorithm
// (c) Cao Siqin 2013.9.26 (first revision)
// (c) Cao Siqin 2013.10.9 (revision) : inversion revised
// (c) Cao Siqin 2014.10.13 (base version)
// (c) Cao Siqin 2018.3.6 (overhaul) : simplified class matrix
// (c) Cao Siqin 2018.3.9 (revision) : MatrixDataContainer::ip auto restore
//          auto restore of ip only occurs at: matrix::(=,+=,-=,*=,/=)
//          auto restore of ip will release: matrix generated in *,det,inv,...
//          auto restore of ip will not prevent the normal allocation of matrix
//          manually recover/restore ip with push_ip/pop_ip
// (c) Cao Siqin 2018.3.12 (revision) : add reference return of element via matrix::a(int,int)
//          This could take the place of the old Matrix::a[][]
// (c) Cao Siqin 2018.3.13 (overhaul) : fix bugs of 2018.3.9
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

namespace MatrixNS {
    const char * Matrix_edition = "(c) Cao Siqin 2013.10.10";
    #ifndef MCONTAINER_SIZE_DEFAULT
        #define MCONTAINER_SIZE_DEFAULT 400000
    #endif
    #ifndef MATRIX_DEFAULT_DIM
        #define MATRIX_DEFAULT_DIM 2
    #endif
    #ifndef STACK_SIZE
        #define STACK_SIZE 10000
    #endif
    #define MIN(a,b) (a<b?a:b)
    #define MAX(a,b) (a>b?a:b)
    class MatrixDataContainer {
      private:
        int ip; int page;
        int sp, stack[STACK_SIZE];
        bool alloc_permanent;
      public:
        char * root;
        int size; int ip_save;
        MatrixDataContainer * next;
      public:
        void init(int size_=MCONTAINER_SIZE_DEFAULT){
            size = size_; //if (size<MCONTAINER_SIZE_DEFAULT) size=MCONTAINER_SIZE_DEFAULT;
            root = 0; next = 0; ip = 0; sp = 0; ip_save = 0; page = 1;
            root = (char*)malloc(sizeof(char)*size);
            alloc_permanent = true;
        }
        void enlarge(int size_=MCONTAINER_SIZE_DEFAULT){
            //if (size_<MCONTAINER_SIZE_DEFAULT) size_=MCONTAINER_SIZE_DEFAULT;
            next = (MatrixDataContainer *)malloc(sizeof(MatrixDataContainer)+sizeof(char)*size_);
            memset(next, 0, sizeof(MatrixDataContainer));
            next->size = size_; next->ip = 0; next->sp = sp;
            next->root = (char*)(((char*)next) + sizeof(MatrixDataContainer));
            next->next = 0; next->page = page + 1;
        }
        void dispose(bool free_root=true){
            if (next){ next->dispose(false); ::free(next); }
            if (free_root) ::free(root);
        }
        char * allocate(int count){
            //fprintf(stderr, "allocate %d from %d\n", count, ip);
            if (ip+count < size){
                  //printf("allocate at %d.%d ", page, ip/8);
                int ir = ip; ip += count; if (alloc_permanent) ip_save = ip;
                  //printf("end at %d.%d\n", page, ip/8);
                for (int i=ir; i<ip; i++) root[i] = 0;
                return &root[ir];
            } else {
                if (!next) enlarge(MAX(size, count));
                return next->allocate(count);
            }
        }
        void reset(){
            sp = 0; ip = 0; ip_save = 0;
            if (next) next->reset(); ip = 0; for (int i=0; i<size; i++) root[i] = 0;
        }
      public:
        void push_ip(){
            if (sp<STACK_SIZE){ stack[sp] = ip; sp ++; }
            if (next) next->push_ip();
        }
        void pop_ip(){
            if (sp>0){ sp --; ip = stack[sp]; }
            if (next) next->pop_ip();
        }
        int getip(){ if (next) return next->getip(); else return ip; }
        int getipsave(){ if (next) return next->getipsave(); else return ip_save; }
        int getpage(){ if (next) return next->getpage(); else return page; }
        void _saveip(){ ip_save = ip; if (next) next->_saveip(); }
        void _restoreip(){ ip = ip_save; if (ip<0) ip = 0; if (next) next->_restoreip(); }
      public:
       // if allocate is temporary, then it will be released at matrix::(=,+=,-=,*=,/=)
        void _set_alloc_permanent(){ alloc_permanent = true; }
        void _set_alloc_temp(){ alloc_permanent = true; }
        bool _is_alloc_permanent(){ return alloc_permanent; }
    };
}

namespace MatrixNS {
    class matrix {
      public:
        int     n;  //dimension
        double* m;  //data array
      protected:
        double* m_new;
        double  __m[MATRIX_DEFAULT_DIM*MATRIX_DEFAULT_DIM];
        MatrixDataContainer * __c; int page, ip;
    public:
        void set_container(MatrixDataContainer * cont){ __c = cont; }
        void init(MatrixDataContainer * cont, int n_){ init(n_, cont); }
        void init(int n, double * m){ __c = 0; this->n = n; this->m = m; m_new = 0; }
        void init(int n, MatrixDataContainer * cont=null){ this->n = n; __c = cont;
            if (n<=MATRIX_DEFAULT_DIM){
                m = __m;
            } else {
                int lend = sizeof(double)*n*n;
                if (__c){ page = __c->getpage(); ip = __c->getip(); }
                char * p = __c?__c->allocate(lend):(char*)malloc(lend);
                m = m_new = (double*) p;
            }
            for (int i=0; i<n*n; i++) m[i] = 0;
        }
        void harness(double * data, int n, MatrixDataContainer * mc=NULL){
            m = data; this->n = n; if (mc) __c = mc;
        }
        void dispose(){ if (!__c&&m_new&&m==m_new) free(m_new); this->n = 0; this->m = this->m_new = NULL; }
        ~matrix(){ dispose(); }
        matrix(){ init(0); }
        matrix(MatrixDataContainer * cont){ init(cont, 1); }
        matrix(MatrixDataContainer * cont, int n_){ init(cont, n_); }
        matrix(int n){ init(n); }
        matrix(int n, double * m){ init(n, m); }
        //matrix(const matrix & o){ memcpy(this, &o, sizeof(matrix)); } // copy constructor ignored because frequently called
      //basic operations
        matrix operator =(double lambda){
            for (int i=0; i<n; i++) for (int j=0; j<n; j++){
                if (i==j) *e(i,j) = lambda; else *e(i,j) = 0;
            }
            return *this;
        }
        matrix operator =(double * array){ for (int i=0; i<n*n; i++) m[i] = array[i]; return *this; }
        matrix operator =(matrix o){
            //printf("(%d.%d,%dx%d).matrix=(%d.%d,%dx%d)%s(ip=%d,save:%d)\n", page,ip/8, n,n, o.page,o.ip/8, o.n,o.n, __c?__c->_is_alloc_permanent()?" P":" T":" T", __c?__c->getip()/8:0, __c?__c->getipsave()/8:0);
            if (__c && __c->_is_alloc_permanent()){
                __c->_restoreip(); int page_c = __c->getpage(); int ip_c = __c->getip(); // free all temporary memory allocated
                if (page>page_c || (page==page_c&&ip>ip_c)){
                    //printf("  matrix.=: %d.%d move to %d.%d", page, ip/8, page_c, ip_c/8);
                    matrix move_to; move_to.init(o.n, __c); memcpy(&move_to, this, sizeof(matrix));
                }
            }

            if (!m || n<o.n){
                if (!o.__c && m_new && m==m_new) dispose();
                init(o.n, o.__c);
            }
            int dim = MIN(n,o.n); for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) *e(i,j) = *o.e(i,j);
            memcpy(__m, o.__m, sizeof(__m));

            return *this;
            
            /*
            if (__c && __c->_is_alloc_permanent()) __c->_restoreip(); // free all temporary memory allocated
            if (!m || o.n<n){ dispose(); memcpy(this, &o, sizeof(matrix)); init(__c, n); }
            int dim = MIN(n,o.n); for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) *e(i,j) = *o.e(i,j);
            return *this;
            */
        }
      //inline methods
        double & a(int a, int b){ return m[a*n+b]; }
        inline double * e(int a, int b){ return &m[a*n+b]; }
      //basic algorithms
        matrix operator +=(matrix o){
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++) *e(i,j) += *o.e(i,j);
            return *this;
        }
        matrix operator -=(matrix o){
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++) *e(i,j) -= *o.e(i,j);
            if (__c && __c->_is_alloc_permanent()) __c->_restoreip(); // free all temporary memory allocated
            return *this;
        }
        matrix operator *=(double l){
            for (int i=0; i<n; i++) for (int j=0; j<n; j++) *e(i,j) *= l;
            return *this;
        }
        matrix operator /=(double l){
            for (int i=0; i<n; i++) for (int j=0; j<n; j++) *e(i,j) /= l;
            return *this;
        }
        matrix operator +(double l){
            if (__c) __c->_set_alloc_temp();
            matrix ret(__c, n);
            for (int i=0; i<n; i++) for (int j=0; j<n; j++){
                *ret.e(i,j) = *e(i,j) + (i==j? l : 0);
            }
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        matrix operator +(matrix o){
            if (__c) __c->_set_alloc_temp();
            matrix ret(__c, MIN(n,o.n));
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *ret.e(i,j) = *e(i,j) + *o.e(i,j);
            }
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        matrix operator -(double l){
            if (__c) __c->_set_alloc_temp();
            matrix ret(__c, n);
            for (int i=0; i<n; i++) for (int j=0; j<n; j++){
                *ret.e(i,j) = *e(i,j) - (i==j? l : 0);
            }
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        matrix operator -(matrix o){
            if (__c) __c->_set_alloc_temp();
            matrix ret(__c, MIN(n,o.n));
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *ret.e(i,j) = *e(i,j) - *o.e(i,j);
            }
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        matrix operator *(double ol){
            if (__c) __c->_set_alloc_temp();
            matrix ret(__c, n);
            for (int i=0; i<n; i++) for (int j=0; j<n; j++){
                *ret.e(i,j) = *e(i,j) * ol;
            }
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        matrix operator /(double ol){
            matrix ret(__c, n);
            for (int i=0; i<n; i++) for (int j=0; j<n; j++){
                *ret.e(i,j) = *e(i,j) / ol;
            }
            return ret;
        }
        matrix product(matrix o, matrix out){
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *out.e(i,j) = 0;
                for (int k=0; k<MIN(n,o.n); k++) *out.e(i,j) += *e(i,k) * *o.e(k,j);
            }
            return out;
        }
        matrix operator *(matrix o){
            if (__c) __c->_set_alloc_temp();
            matrix ret(__c, MIN(n,o.n));
            for (int i=0; i<MIN(n,o.n); i++) for (int j=0; j<MIN(n,o.n); j++){
                *ret.e(i,j) = 0;
                for (int k=0; k<MIN(n,o.n); k++) *ret.e(i,j) += *e(i,k) * *o.e(k,j);
            }
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        matrix operator /(matrix ol){
            return *this * ol.inverse();
        }
      //Transpose
        matrix transpose(){
            for (int i=0; i<n; i++) for (int j=i; j<n; j++){
                double t = *e(i,j); *e(i,j) = *e(j,i); *e(j,i) = t;
            }
            return * this;
        }
        matrix tr(){ return transpose(); return * this; }
      //Trace
        double trace(){
            double tr = 0; for (int i=0; i<n; i++) tr += *e(i,i); return tr;
        }
      public:
       //line transformation: basic transformation
        void row_weight_add_to(int a, int b, double k, int c0=0, int cn=-1){ //row a product k then add to row b (B += k*A), from col c0 to cn
            if (cn<c0) cn = n-1; for (int i=c0; i<=cn; i++) *e(b, i) += *e(a, i) * k;
        }
        void col_weight_add_to(int a, int b, double k, int c0=0, int cn=-1){ //col a product k then add to col b (B += k*A), from row c0 to cn
            if (cn<c0) cn = n-1; for (int i=c0; i<=cn; i++) *e(i, b) += *e(i, a) * k;
        }
      public:
       //trianglization
        bool renom_diagnol_row(matrix oi){
            /*for (int k=0; k<n; k++) if (*e(k,k)==0){
                int l=-1; for (int kk=0; kk<n; kk++){ if (*e(kk,k)!=0){ l = kk; break; } }
                if (l<0) return false;
            }*/
            for (int k=0; k<n; k++) if (*e(k,k)==0){
                int l=-1; for (int kk=0; kk<n; kk++){ if (*e(kk,k)!=0){ l = kk; break; } }
                if (l<0) return false;
                row_weight_add_to(l, k, 1);
                oi.row_weight_add_to(l, k, 1);
            }
            return true;
        }
        matrix up_triangle(matrix oi){ // upper trianglize with row operations (oi^-1 * up_triangle = original)
            bool reversible = true;
            reversible = renom_diagnol_row(oi);
            if (reversible) for (int u=0; u<=n-2; u++){
                for (int k=u+1; k<=n-1; k++){
                  if (*e(u,u) == 0){
                    reversible = false;
                  } else {
                    double a = - *e(k,u) / *e(u,u);
                    row_weight_add_to(u, k, a, u+1); *e(k,u) = 0;
                    oi.row_weight_add_to(u, k, a);
                  }
                }
            }
            if (!reversible) return matrix(0);
            return *this;
        }
        
        bool renom_diagnol_col(matrix oi){
            /*for (int k=0; k<n; k++) if (*e(k,k)==0){
                int l=-1; for (int kk=0; kk<n; kk++){ if (*e(kk,k)!=0){ l = kk; break; } }
                if (l<0) return false;
            }*/
            for (int k=0; k<n; k++) if (*e(k,k)==0){
                int l=-1; for (int kk=0; kk<n; kk++){ if (*e(k,kk)!=0){ l = kk; break; } }
                if (l<0) return false;
                col_weight_add_to(l, k, 1);
                oi.col_weight_add_to(l, k, 1);
            }
            return true;
        }
        matrix low_triangle(matrix oi){ // lower trianglize with column operations (down_triangle / oi = original)
            bool reversible = true;
            reversible = renom_diagnol_col(oi);
            if (reversible) for (int u=0; u<=n-2; u++){
                for (int k=u+1; k<=n-1; k++){
                  if (*e(u,u) == 0){
                    reversible = false;
                  } else {
                    double a = - *e(u,k) / *e(u,u);
                    col_weight_add_to(u, k, a, u+1); *e(u,k) = 0;
                    oi.col_weight_add_to(u, k, a);
                  }
                }
            }
            if (!reversible) return matrix(0);
            return *this;
        }
       //diagnol renormalization: self renormalization, eleminate diagnol 0 element
        double determin(bool * reversible = null){
            if (__c) __c->_set_alloc_temp();
            matrix tmp(__c, n); tmp = *this;
            matrix tmp2 = tmp.up_triangle(matrix(0));
              //if (tmp2.n!=tmp.n){ if (reversible) *reversible = false; else *reversible = true; } //ƒÊ≤ª¥Ê‘⁄
            double ret = 1; for (int i=0; i<n; i++) ret *= *tmp2.e(i,i);
            if (__c) __c->_set_alloc_permanent();
            return ret;
        }
        double det(bool * reversible = null){ return determin(reversible); }
        double determinant(bool * reversible = null){ return determin(reversible); }
      //Inversion
        matrix inverse(){
            if (__c) __c->_set_alloc_temp();
            matrix oi(__c, n);  matrix ol(__c, n);
            if (!inverse(oi, ol)) oi = 0.0;
            if (__c) __c->_set_alloc_permanent();
            return oi;
        }
        bool inverse(matrix & oi, matrix & ol){ //ol is a temporary matrix, which can be *this
            if (__c) __c->_set_alloc_temp();
            ol = *this; oi = 1;
          //1&2 lower triangle transformation: lower triangle transformed to 0
            matrix t = ol.up_triangle(oi); if (t.n!=ol.n) return false;  //no inversion
          //3. check the diagnol
            //for (int i=0; i<ol.n; i++)  if (*ol.e(i,i)==0) return false;  //no inversion
          //4. upper triangle transformation: upper triangle transformed to 0
            for (int u=ol.n-1; u>=1; u--){
                for (int k=u-1; k>=0; k--){
                    double a = - *ol.e(k,u) / *ol.e(u,u);
                    ol.row_weight_add_to(u, k, a, 0, u-1); *ol.e(k,u) = 0;
                    oi.row_weight_add_to(u, k, a);
                }
            }
          //5. renormalized with diagnol elements
            for (int l=0; l<ol.n; l++){
                double factor = *ol.e(l,l);
                for (int i=0; i<ol.n; i++){
                    *ol.e(l,i) /= factor;
                    *oi.e(l,i) /= factor;
                }
            }
            if (__c) __c->_set_alloc_permanent();
            return true;
        }
    };
}

