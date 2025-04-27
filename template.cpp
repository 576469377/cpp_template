#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

typedef pair<int, int> PII;
typedef long long LL;

const int INF = 0x3f3f3f3f;
const int N = 100010;

/////////////////////////////////////////////////////////

// 快速排序
void quick_sort(int q[], int l, int r) {
    if (l >= r) return;
    int x = q[(l + r) >> 1];
    int i = l - 1, j = r + 1;
    while (i < j) {
        while (q[++i] < x);
        while (q[--j] > x);
        if (i < j) swap(q[i], q[j]);
    }
    quick_sort(q, l, j), quick_sort(q, j + 1, r);
}

/////////////////////////////////////////////////////////

// 归并排序
void merge_sort(int q[], int tmp[], int l, int r) {
    if (l >= r) return;
    int mid = (l + r) >> 1;
    merge_sort(q, tmp, l, mid), merge_sort(q, tmp, mid + 1, r);
    int i = l, j = mid + 1, k = 0;
    while (i <= mid && j <= r) {
        if (q[j] < q[i]) tmp[k++] = q[j++];
        else tmp[k++] = q[i++];
    }
    while (i <= mid) tmp[k++] = q[i++];
    while (j <= r) tmp[k++] = q[j++];
    for (i = l, j = 0; i <= r; ++i, ++j) q[i] = tmp[j];
}

/////////////////////////////////////////////////////////

// 整数二分
// 区间[l, r]被划分成[l, mid]和[mid + 1, r]
bool check(int x) {}
int b_search_1(int l, int r) {
    while (l < r) {
        int mid = (l + r) >> 1;
        if (check(mid)) r = mid;
        else l = mid + 1;
    }
    return l;
}

// 区间[l, r]被划分成[l, mid - 1]和[mid, r]
int b_search_2(int l, int r) {
    while (l < r) {
        int mid = (l + r + 1) >> 1;
        if (check(mid)) l = mid;
        else r = mid - 1;
    }
    return l;
}

// 浮点数二分
bool check(double x) {}
const double eps = 1e-8;
double b_search_3(double l, double r) {
    while (r - l > eps) {
        double mid = (l + r) / 2;
        if (check(mid)) r = mid;
        else l = mid;
    }
    return l;
}

/////////////////////////////////////////////////////////

// string数字转换为vector倒序数字
vector<int> string_to_vector(string a) {
    vector<int> ans;
    for (int i = a.size() - 1; i >= 0; --i) ans.push_back(a[i] - '0');
    return ans;
}

// 高精度加法
// a, b均为倒序
vector<int> add(vector<int> &a, vector<int> &b) {
    if (a.size() < b.size()) return add(b, a);
    vector<int> c;
    int t = 0;
    for (int i = 0; i < a.size(); ++i) {
        t += a[i];
        if (i < b.size()) t += b[i];
        c.push_back(t % 10);
        t /= 10;
    }
    if (t > 0) c.push_back(t);
    return c;
}

// 压位：string存为vector
const int BASE = 1e9; // 每9位作为一组
vector<int> string_to_vector_compress(string &a) {
    vector<int> ans;
    for (int i = a.size() - 1, j = 0, s = 0, t = 1; i >= 0; --i) {
        s += (a[i] - '0') * t;
        ++j, t *= 10;
        if (j == 9 || i == 0) { // 满9位或者当前为最后一位，进位
            ans.push_back(s);
            j = 0, s = 0, t = 1;
        }
    }
    return ans;
}

// 打印压位数字
void print_compress(vector<int> &a) {
    printf("%d", a.back());
    for (int i = a.size() - 2; i >= 0; --i) printf("%09d", a[i]);
    cout << endl;
}

// 高精度加法（压位）
vector<int> add_compress(vector<int> &a, vector<int> &b) {
    if (a.size() < b.size()) return add_compress(b, a);
    vector<int> c;
    int t = 0;
    for (int i = 0; i < a.size(); ++i) {
        t += a[i];
        if (i < b.size()) t += b[i];
        c.push_back(t % BASE);
        t /= BASE;
    }
    if (t > 0) c.push_back(t);
    return c;
}

// 比较两个vector数字是否a >= b
bool cmp_vector_vector(vector<int> &a, vector<int> &b) {
    if (a.size() != b.size()) return a.size() > b.size();
    for (int i = a.size() - 1; i >= 0; --i) {
        if (a[i] != b[i]) return a[i] > b[i];
    }
    return true;
}

// 高精度减法，默认a >= b
vector<int> sub(vector<int> &a, vector<int> &b) {
    vector<int> c;
    int t = 0;
    for (int i = 0; i < a.size(); ++i) {
        t = a[i] - t;
        if (i < b.size()) t -= b[i];
        if (t >= 0) {
            c.push_back(t);
            t = 0;
        } else {
            c.push_back(t + 10);
            t = 1;
        }
    }
    while (c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}

// 高精度乘法
vector<int> mul(vector<int> &a, int b) {
    int t = 0;
    vector<int> c;
    for (int i = 0; i < a.size(); ++i) {
        t += a[i] * b;
        c.push_back(t % 10);
        t /= 10;
    }
    while (t) {
        c.push_back(t % 10);
        t /= 10;
    }
    while (c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}


// 高精度除法（a为正序数字，r为余数）
vector<int> div(vector<int> &a, int b, int &r) {
    r = 0;
    vector<int> c;
    for (int i = 0; i < a.size(); ++i) {
        r = r * 10 + a[i];
        c.push_back(r / b);
        r %= b;
    }
    reverse(c.begin(), c.end());
    while (c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}

// 数组高精度加法（其余类似）
void add(LL a[], LL b[]) {
    static LL c[M];
    memset(c, 0, sizeof c);
    for (int i = 0, t = 0; i < M; ++i) {
        t += a[i] + b[i];
        c[i] = t % 10;
        t /= 10;
    }
    memcpy(a, c, sizeof c);
}

/////////////////////////////////////////////////////////

// 2D前缀和
int a[N], s[N];
// 初始化
for (int i = 1; i <= n; ++i) s[i] = s[i - 1] + a[i];
//查询[l, r]的和
int get(int l, int r) {
    return s[r] - s[l - 1];
}

// 3D前缀和
int a[N][N], s[N][N];
// 初始化
for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
        s[i][j] = a[i][j] + s[i][j - 1] + s[i - 1][j] - s[i - 1][j - 1];
    }
}
// 查询矩阵[x1, y1], [x2, y2]的和
int get(int x1, int y1, int x2, int y2) {
    return s[x2][y2] - s[x2][y1 - 1] - s[x1 - 1][y2] + s[x1 - 1][y1 - 1];
}

/////////////////////////////////////////////////////////

// 2D差分
int a[N], b[N];
// [l, r]区间上的所有元素加上c
void add(int l, int r, int c) {
    b[l] += c;
    b[r + 1] -= c;
}
// 初始化
for (int i = 1; i <= n; ++i) add(i, i, a[i]);
// 计算原数组
for (int i = 1; i <= n; ++i) b[i] += b[i - 1];

// 3D差分
int a[N][N], b[N][N];
// [x1, y1], [x2, y2]矩形上所有元素加上c
void add(int x1, int y1, int x2, int y2, int c) {
    b[x1][y1] += c;
    b[x1][y2 + 1] -= c;
    b[x2 + 1][y1] -= c;
    b[x2 + 1][y2 + 1] += c;
}
// 初始化
for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
        add(i, j, i, j, a[i][j]);
    }
}
// 计算原数组
for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
        b[i][j] += b[i - 1][j] + b[i][j - 1] - b[i - 1][j - 1];
    }
}

/////////////////////////////////////////////////////////

// 位运算
// 获取最低一位的1，e.g.101010 -> 10
int lowbit(int x) {
    return x & -x;
}

/////////////////////////////////////////////////////////

// 离散化
// a: 差分数组
// s: 原数组
int a[N], s[N];
vector<int> alls; // 出现过的元素放入alls

// 排序并唯一化
sort(alls.begin(), alls.end());
alls.erase(unique(alls), alls.end());

// 返回>=x的第一个数的位置+1(差分数组下标从1开始)
int find(int x) {
    int l = 0, r = alls.size() - 1;
    while (l < r) {
        int mid = (l + r) >> 1;
        if (alls[mid] >= x) r = mid;
        else l = mid + 1;
    }
    return l + 1;
}

/////////////////////////////////////////////////////////

// 合并区间
void merge(vector<PII> &segs) {
    vector<PII> ans;
    sort(segs.begin(), segs.end());
    int l = -INF, r = -INF;
    for (auto &it: segs) {
        if (r < it.first) {
            if (l != -INF) ans.push_back({l, r});
            l = it.first, r = it.second;
        } else r = max(r, it.second);
    }
    if (l != -INF) ans.push_back({l, r});
    segs = ans;
}

/////////////////////////////////////////////////////////

// 单链表
int head = -1, e[N], ne[N], idx;

void add_to_head(int x) {
    e[idx] = x, ne[idx] = head, head = idx++;
}

void add(int k, int x) { // 在第k个插入的数后面插入一个数（从0开始）
    e[idx] = x, ne[idx] = ne[k], ne[k] = idx++;
}

void remove(int k) { // 删除第k个插入的数后面的数
    ne[k] = ne[ne[k]];
}

/////////////////////////////////////////////////////////

// 双链表
int e[N], l[N], r[N], idx;
// 初始化：0号和1号点作为头尾节点
r[0] = 1, l[1] = 0, idx = 2;

// 第0个点和第1个点已经被使用，在操作insert时候应该要看情况加2
void insert(int k, int x) { // 在第k个插入的数后面插入一个数
    e[idx] = x;
    l[idx] = k, r[idx] = r[k];
    l[r[idx]] = idx, r[k] = idx++;
}

void remove(int k) { // 删除第k个数
    l[r[k]] = l[k];
    r[l[k]] = r[k];
}

/////////////////////////////////////////////////////////

// 栈
int stk[N], tt;

stk[++tt] = x; // 入栈

--tt; // 出栈

stk[tt]; // 栈顶

// 判空 if (tt): 非空

/////////////////////////////////////////////////////////

// 队列
int q[N], hh, tt = -1;

q[++tt] = x; // 入队

++hh; // 出队

// 判空 if (hh <= tt): 非空

/////////////////////////////////////////////////////////

// 单调栈
// 用于计算上一个与当前元素有某种大小关系的元素
// e.g. 计算上一个比当前元素小的元素
//      单调栈的元素单调性满足x_1 < x_2 < ... < x_tt

// 此处栈内存元素
while (tt && stk[tt] >= x) --tt; // 不满足，一直出栈
stk[++tt] = x; // 直到满足，入栈

/////////////////////////////////////////////////////////

// 单调队列
// 用于计算滑动窗口内的某种性质
// e.g. 计算滑动窗口最大值
//      单调队列中的元素单调性满足x_hh > ... > x_tt

// 此处队列内存下标
// 当前元素为a[i]
while (hh <= tt && i - q[hh] >= k) ++hh; // 超出窗口，一直出队
while (hh <= tt && a[q[tt]] <= a[i]) --tt; // 不满足，一直出队
q[++tt] = i;
if (i + 1 >= k) // 获取滑动窗口最大值

/////////////////////////////////////////////////////////

// KMP
int n;
char s[N];
int m;
char p[N];
int ne[N];

// 版本1：下标从0开始
// j表示当前失配后下一个应该匹配的字符下标
void get1() {
    int i = 0, j = ne[0] = -1;
    while (i < n) {
        if (j < 0 || p[i] == p[j]) ++i, ++j, ne[i] = j;
        else j = ne[j];
    }
}
void get1() { // 优化版
    int i = 0, j = ne[0] = -1;
    while (i < n) {
        if (j < 0 || p[i] == p[j]) ++i, ++j, ne[i] = (p[i] == p[j] ? ne[j] : j);
        else j = ne[j];
    }
}
// 匹配
int i = 0, j = 0;
while (i < m && j < n) {
    if (j < 0 || s[i] == p[j]) ++i, ++j;
    else j = ne[j];
    if (j == n) {
        cout << i - j << " ";
        j = ne[j];
    }
}

// 版本2：下标从1开始（计算的是前缀函数）
// j表示当前成功匹配的长度
void get2() { // 字符串下标从1开始
    for (int i = 2, j = 0; i <= n; ++i) {
        while (j && p[i] != p[j + 1]) j = ne[j];
        if (p[i] == p[j + 1]) ++j;
        ne[i] = j;
    }
}
// 匹配
for (int i = 1, j = 0; i <= m; ++i) {
    while (j && s[i] != p[j + 1]) j = ne[j];
    if (s[i] == p[j + 1]) ++j;
    if (j == n) {
        cout << i - j << " ";
        j = ne[j];
    }
}   

/////////////////////////////////////////////////////////

// Trie

// 空节点和根节点的值都是0
int son[N][26], cnt[N], idx;
char str[N];

void insert(char *str) {
    int p = 0;
    for (int i = 0; str[i]; ++i) {
        int u = str[i] - 'a';
        if (!son[p][u]) son[p][u] = idx++;
        p = son[p][u];
    }
    ++cnt[p];
}

// 查找str出现的次数
int query(char *str) {
    int p = 0;
    for (int i = 0; str[i]; ++i) {
        int u = str[i] - 'a';
        if (!son[p][u]) return 0;
        p = son[p][u];
    }
    return cnt[p];
}

/////////////////////////////////////////////////////////

// 并查集
// 如果需要可以维护一个sz[N]作为连通块大小
int p[N];

// 初始化
for (int i = 1; i <= n; ++i) p[i] = i;

// 路径压缩
int find(int x) {
    return x == p[x] ? x : p[x] = find(p[x]);
}

// 合并a, b
p[fnid(a)] = find(b);

/////////////////////////////////////////////////////////

// 堆
// 下滤
void down(int u) {
    int t = u;
    if (u * 2 <= n && h[u * 2] < h[t]) t = u * 2;
    if (u * 2 + 1 <= n && h[u * 2 + 1] < h[t]) t = u * 2 + 1;
    if (u != t) swap(h[u], h[t]), down(t);
}

// 上滤
void up(int u) {
    while (u / 2 && h[u / 2] > h[u]) swap(h[u / 2], h[u]), u /= 2;
}

// 建堆（O(n)）
for (int i = n / 2; i; --i) down(i);

// 插入
h[++n] = x;
up(n);

// 删除
h[1] = h[n--];
down(1);

/////////////////////////////////////////////////////////

// 支持对第k个插入的元素进行操作的堆

// ph[k] = i: 第k个插入的元素对应的堆的下标为i
// hp[i] = k: 下标为i的堆中元素是第k个插入的
// n: 当前堆中的元素个数
// cnt: 当前即将插入的元素是第几个

int n, cnt;
int h[N];
int ph[N], hp[N];

// 交换下标为a, b的两个元素
void h_swap(int a, int b) {
    swap(ph[hp[a]], ph[hp[b]]);
    swap(hp[a], hp[b]);
    swap(h[a], h[b]);
}

void down(int u) {
    int t = u;
    if (u * 2 <= n && h[u * 2] < h[t]) t = u * 2;
    if (u * 2 + 1 <= n && h[u * 2 + 1] < h[t]) t = u * 2 + 1;
    if (u != t) h_swap(u, t), down(t);
}

void up(int u) {
    while (u / 2 && h[u / 2] > h[u]) h_swap(u / 2, u), u /= 2;
}

// 插入
++n, ++cnt;
ph[cnt] = n, ph[n] = cnt;
h[n] = x;
up(n);

// 删除堆顶
h_swap(1, n--);
down(1);

// 删除第k个插入的元素
int t = ph[k];
h_swap(t, n--);
up(t), down(t);

// 修改第k个插入的元素
int t = ph[k];
h[t] = x;
up(t), down(t);

/////////////////////////////////////////////////////////

// 哈希表
// 开放寻址（闭散列）

// 开放寻址表长一般取作数据范围的2倍以上的素数（此种方式寻址只能使用到一半左右的桶）
// 使用一个很大的数字作为空标记（null = 0x3f3f3f3f）

const int null = 0x3f3f3f3f;

int h[N];

// 查找（失败返回null所在的位置）
int find(int x) {
    int t = (x % N + N) % N;
    while (h[t] != null && h[t] != x) {
        ++t; // 线性试探
        if (t == N) t = 0; // 超出则重新回到开头
    }
    return t;
}

// 初始化
// memset按字节赋值，每一个字节都赋上3f，相当于每一个数都赋为0x3f3f3f3f
memset(h, 0x3f, sizeof h);

// 插入
h[find(x)] = x;

// 查找成功 if (h[find(x)] != null)

/////////////////////////////////////////////////////////

// 哈希表
// 封闭寻址（拉链法，开散列）

// 表长取作大于数据范围的第一个素数
const int N = 100003;

int h[N], e[N], ne[N], idx;

// 插入
void insert(int x) {
    int k = (x % N + N) % N;
    e[idx] x, ne[idx] = h[k], h[k] = idx++;
}

bool find(int x) {
    int k = (x % N + N) % N;
    for (int i = h[k]; ~i; i = ne[i])
        if (e[i] == x) return true;
    return false;
}

/////////////////////////////////////////////////////////

// 字符串哈希

typedef unsigned long long ULL;

// 将字符串转换为P进制数（131为经验值）
const int P = 131;

char str[N];
ULL h[N], p[N]; // h为前缀哈希值，p为预处理的当前位的权重

// 得到str[l:r]的哈希值
ULL get(int l, int r) {
    return h[r] - h[l - 1] * p[r - l + 1]; // 使用ULL（模2^64，一旦溢出相当于自动取模）
}

// 预处理h和p
p[0] = 1;
for (int i = 1; i <= n; ++i) {
    h[i] = h[i - 1] * P + str[i];
    p[i] = p[i - 1] * P;
}

/////////////////////////////////////////////////////////

// bfs搜索方格图的最短路
// 从 (0, 0) 走到 (n - 1, m - 1) 的最短路（移动一次的距离为 1）

const int dirs[4][2] = {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};

int n, m;
int dis[N][N];

int bfs() {
    queue<int> q;
    memset(dis, -1, sizeof dis);
    dis[0][0] = 0;
    q.push(0);
    while (q.size()) {
        int f = q.front(); q.pop();
        int x = f / m, y = f % m;
        for (int i = 0; i < 4; ++i) {
            int xx = x + dirs[i][0], yy = y + dirs[i][1];
            if (xx >= 0 && xx < n && yy >= 0 && yy < m && dis[xx][yy] == -1) {
                dis[xx][yy] = dis[x][y] + 1;
                q.push(xx * m + yy);
            }
        }
    }
    return dis[n - 1][m - 1];
}

/////////////////////////////////////////////////////////

// bfs搜索图的最短路

int h[N], e[N], ne[N], idx;
int dis[N];

void add(int a, int b) {
    e[idx] = b, ne[idx] = h[a], h[a] = idx++;
}

int bfs() {
    dis[1] = 0;
    queue<int> q;
    q.push(1);
    while (q.size()) {
        int f = q.front(); q.pop();
        if (f == n) return dis[f];
        for (int i = h[f]; ~i; i = ne[i]) {
            int j = e[i];
            if (dis[j] == -1) {
                dis[j] = dis[f] + 1;
                q.push(j);
            }
        }
    }
    return -1;
}

/////////////////////////////////////////////////////////

// 0入度拓扑排序

int n, m;
int h[N], e[N], ne[N], idx;
int in[N];

void add(int a, int b) {
    e[idx] = b, ne[idx] = h[a], h[a] = idx++;
    ++in[b];
}

void top_sort() {
    vector<int> res;
    queue<int> q;
    for (int i = 1; i <= n; ++i) if (in[i] == 0)
        q.push(i);
    while (q.size()) {
        int f = q.front(); q.pop();
        res.push_back(f);
        for (int i = h[f]; ~i; i = ne[i]) {
            int j = e[i];
            if (--in[j] == 0) {
                q.push(j);
            }
        }
    }
    if (res.size() == n)
        for (int i = 0; i < n; ++i) cout << res[i] << " ";
    else cout << "-1";
    cout << endl;
}

/////////////////////////////////////////////////////////

// 朴素dijkstra
// 适用于稠密图
// 使用邻接矩阵

int n, m;
int g[N][N], dis[N], st[N];

// 初始化
memset(g, 0x3f, sizeof g);
g[x][y] = min(g[x][y], z); // 重边算短的（读入z）

int dijkstra() {
    memset(dis, 0x3f, sizeof dis);
    dis[1] = 0;
    for (int i = 0; i < n; ++i) {
        // 找到当前已经确定的最近的点
        int u = -1;
        for (int j = 1; j <= n; ++j)
            if (!st[j] && (u < 0 || dis[j] < dis[u]))
                u = j;
        st[u] = 1;
        // 更新下一个点（更新从当前点到其余所有点各自的距离）
        for (int v = 1; v <= n; ++v)
            dis[v] = min(dis[v], dis[u] + g[u][v]);
    }
    if (dis[n] == 0x3f3f3f3f) return -1;
    return dis[n];
}

/////////////////////////////////////////////////////////

// 堆优化dijkstra
// 适用于稀疏图
// 使用邻接表

typedef pair<int, int> PII; // dis, pos

int n, m;
int h[N], e[N], ne[N], idx;
int w[N];
int dis[N], st[N];

void add(int a, int b, int c) {
    e[idx] = b, w[idx] = c, ne[idx] = h[a], h[a] = idx++;
}

int dijkstra() {
    memset(dis, 0x3f, sizeof dis);
    dis[1] = 0;
    priority_queue<PII, vector<PII>, greater<PII>> pq;
    pq.push({0, 1});
    while (pq.size()) {
        int u = pq.top().second; pq.pop();
        if (st[u]) continue;
        st[u] = 1;
        for (int i = h[u]; ~i; i = ne[i]) {
            int v = e[i];
            if (dis[v] > dis[u] + w[i]) {
                dis[v] = dis[u] + w[i];
                pq.push({dis[v], v});
            }
        }
    }
    if (dis[n] == 0x3f3f3f3f) return -1;
    return dis[n];
}

/////////////////////////////////////////////////////////

// Bellman-Ford
// 适用于存在负权边求最短路，可以给出源点到终点的边数限制
// 时间复杂度： O(mn)

// N点数，M边数
const int N = 510, M = 10010;

struct E {
    int a, b, c;
}es[M];

int n, m, k;
int dis[N], last[N];

int bellman_ford() {
    memset(dis, 0x3f, sizeof dis);
    dis[1] = 0;
    for (int i = 0; i < k; ++i) {
        memcpy(last, dis, sizeof dis);
        for (int j = 0; j < m; ++j) {
            E e = es[j];
            dis[e.b] = min(dis[e.b], last[e.a] + e.c);
        }
    }
    // 当存在负权边时，无穷大距离可能会被更新为比无穷大小的一个数
    if (dis[n] > 0x3f3f3f3f - N * M) return 0x3f3f3f3f;
    return dis[n];
}

/////////////////////////////////////////////////////////

// spfa
// 使用队列对Bellman-Ford算法进行优化即得到spfa算法（只有被更新的点才可能入队）
// 时间复杂度：平均O(m)，最差O(mn)

int n, m;
int h[N], e[N], ne[N], idx;
int w[N];
int dis[N], st[N]; // st表示一个元素在不在队列里

void add(int a, int b, int c) {
    e[idx] = b, w[idx] = c, ne[idx] = h[a], h[a] = idx++;
}

int spfa() {
    memset(dis, 0x3f, sizeof dis);
    queue<int> q;
    q.push(1);
    dis[1] = 0;
    st[1] = 1;
    while (q.size()) {
        int u = q.front(); q.pop();
        st[u] = 0; // 出队置0
        for (int i = h[u]; ~i; i = ne[i]) {
            int v = e[i];
            if (dis[v] > dis[u] + w[i]) {
                dis[v] = dis[u] + w[i];
                if (!st[v]) { // 被更新的点才可能入队
                    st[v] = 1;
                    q.push(v);
                }
            }
        }
    }
    if (dis[n] > 0x3f3f3f3f / 2) return 0x3f3f3f3f;
    return dis[n];
}

/////////////////////////////////////////////////////////

// spfa判断负环
// 如果到达当前点的边大于等于n，说明存在负权环路
// 1. 边数等于n，则点数为n+1，存在环路
// 2. 如果环路是正权，则不可能更新到n条边（总权值减少的环路才可能反复更新，负权环路满足了这一性质）

int n, m;
int h[N], e[M], ne[M], idx;
int w[M];
int dis[N], st[N];
int cnt[N]; // 到达当前点的最短路的边数

int spfa() {
    memset(dis, 0x3f, sizeof dis);
    queue<int> q;
    // 将所有点放入队列，环路可能以任何一个点为起点
    for (int i = 1; i <= n; ++i) {
        st[i] = 1;
        q.push(i);
    }
    while (q.size()) {
        int u = q.front(); q.pop();
        st[u] = 0;
        for (int i = h[u]; ~i; i = ne[i]) {
            int v = e[i];
            if (dis[v] > dis[u] + w[i]) {
                dis[v] = dis[u] + w[i];
                cnt[v] = cnt[u] + 1;
                if (cnt[v] >= n) return 1;
                if (!st[v]) {
                    st[v] = 1;
                    q.push(v);
                }
            }
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////

// Floyd算法
// 用于求解多源最短路
// 时间复杂度：O(n^3)

int dis[N][N];

void floyd() {
    for (int k = 1; k <= n; ++k) {
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                dis[i][j] = min(dis[i][j], dis[i][k] + dis[k][j]);
            }
        }
    }
}

// 初始化
for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= n; ++j) {
        if (i != j) dis[i][j] = 0x3f3f3f3f;
        else dis[i][j] = 0; // 自环
    }
}

/////////////////////////////////////////////////////////

// 朴素Prim算法
// 用于稠密图求MST
// 时间复杂度：O(n^2)

int n, m;
int g[N][N];
int dis[N], st[N];

// MST权重之和
int prim() {
    int ans = 0;
    memset(dis, 0x3f, sizeof dis);
    for (int i = 0; i < n; ++i) {
        int u = -1;
        for (int j = 1; j <= n; ++j) {
            if (!st[j] && (u < 0 || dis[j] < dis[u])) u = j;
        }
        st[u] = 1;
        // i == 1时第一条边加入
        if (i && dis[u] == INF) return INF;
        if (i) ans += dis[u];
        for (int v = 1; v <= n; ++v) {
            dis[v] = min(dis[v], g[u][v]);
        }
    }
    return ans;
}

// 初始化（无向图）
memset(g, 0x3f, sizeof g);
g[a][b] = g[b][a] = min(g[a][b], c);

/////////////////////////////////////////////////////////

// Kruskal算法
// 用于稀疏图求MST
// 时间复杂度：O(mlogm)

int n, m;
int p[N];

struct E {
    int a, b, w;
    bool operator<(const E &rhs) const {
        return w < rhs.w;
    }
}es[M];

int find(int x) {
    return x == p[x] ? x : p[x] = find(p[x]);
}

int kruskal() {
    int ans = 0, cnt = 0; // cnt记录边的数量
    for (int i = 1; i <= n; ++i) p[i] = i;
    sort(es, es + m);
    for (int i = 0; i < m; ++i) {
        int a = find(es[i].a), b = find(es[i].b);
        if (a != b) {
            p[a] = b;
            ans += es[i].w;
            ++cnt;
        }
    }
    if (cnt == n - 1) return ans;
    return -1;
}

/////////////////////////////////////////////////////////

// 判断二分图（染色法）
// dfs染色
// 一个图是二分图的充要条件是该图不存在奇数环
// 时间复杂度：O(m + n)

int n, m;
int h[N], e[M], ne[M], idx;
int color[N];

void add(int a, int b) {
    e[idx] = b, ne[idx] = h[a], h[a] = idx++;
}

int dfs(int u, int c) {
    color[u] = c;
    for (int i = h[u]; ~i; i = ne[i]) {
        int v = e[i];
        if (!color[v] && !dfs(v, 3 - c)) return 0;
        else if (color[v] == c) return 0;
    }
    return 1;
}

int is_bg() {
    int ans = 1;
    for (int i = 1; i <= n; ++i)
        if (!color[i] && !dfs(i, 1)) { // 对于每一个连通块均从1开始染色
            ans = 0;
            break;
        }
    return ans;
}

add(a, b), add(b, a); // 无向图初始化

/////////////////////////////////////////////////////////

// 匈牙利算法
// 多项式时间之内寻找最优化任务分配（二分图匹配）
// 时间复杂度：O(mn)，平均远小于O(mn)

int n1, n2, m;
int h[N], e[M], ne[M], idx;
int match[N], st[N];

int find(int u) {
    for (int i = h[u]; ~i; i = ne[i]) {
        int v = e[i];
        if (!st[v]) {
            st[v] = 1;
            // 没匹配，则直接匹配
            // 匹配，则让匹配的的点去匹配下一个（match[v]把v让给u）
            if (!match[v] || find(match[v])) {
                match[v] = u;
                return 1;
            }
        }
    }
    return 0;
}

int bg_match() {
    int ans = 0;
    for (int i = 1; i <= n1; ++i) {
        memset(st, 0, sizeof st); // 可能会涉及之前的点让给当前点，因此重置st
        ans += find(i);
    }
    return ans;
}

/////////////////////////////////////////////////////////

// 数学

/////////////////////////////////////////////////////////

// 快速幂
// 时间复杂度：O(logn)，n为指数的大小

LL pow(int a, int b, int p) {
    LL ans = 1 % p;
    while (b) {
        if (b & 1) ans = ans * a % p;
        a = (LL)a * a % p;
        b >>= 1;
    }
    return ans;
}

// 逆元
/*
    欧拉定理：若a与n互质，则a ^ phi(n) = 1(mod n)
    证明：
    令a_1, a_2, ..., a_phi(n)为[1, n]里所有与n互质的数
    由于a与n互质，因此a * a_1, a * a_2, ... , a * a_phi(n)同样为phi(n)个与n互质的数
    即a ^ phi(n) * a_1 * ... * a_phi(n) = a_1 * ... * a_phi(n)(mod n)
    因此a ^ phi(n) = 1(mod n)
    当n为质数时，可以导出费马小定理
    a ^ phi(p) = 1(mod p)
    即a ^ (p - 1) = 1(mod p)

    a存在乘法逆元x的充要条件是a与p互质，ax = 1(mod p)
    由费马小定理可以得到a * a ^ (p - 2) = 1(mod p)
    乘法逆元为a ^ (p - 2)
*/
pow(a, p - 2, p);

/////////////////////////////////////////////////////////

// 试除法判断质数
// 时间复杂度：O(sprt(n))

int is_prime(int x) {
    if (x < 2) return 0;
    for (int i = 2; i <= x / i; ++i) {
        if (x % i == 0) return 0;
    }
    return 1;
}

/////////////////////////////////////////////////////////

// 试除法分解质因数
// 1. 对于当前质因子，连续从x中除去，直到除尽
// 2. 大于sqrt(x)的质因子只可能有一个（如果有两个，则相乘将大于x）

void factorize(int x) {
    for (int i = 2; i <= x / i; ++i) {
        if (x % i == 0) {
            int c = 0;
            while (x % i == 0) {
                ++c;
                x /= i;
            }
            cout << i << " " << c << endl;
        }
    }
    if (x > 1) cout << x << " " << 1 << endl;
    cout << endl;
}

/////////////////////////////////////////////////////////

// 筛法

// primes：存放质数
// cnt：当前的质数个数
// st：当前是否为质数，默认0为质数
int primes[N], cnt;
int st[N];

// 朴素筛（埃氏筛）：时间复杂度O(nloglogn)
void get_primes1(int n) {
    for (int i = 2; i <= n; ++i) {
        if (st[i]) continue; // 合数（已经被筛过）
        primes[cnt++] = i;
        for (int j = i + i; j <= n; j += i) st[j] = 1;
    }
}

// 线性筛（欧拉筛）：时间复杂度O(n)
// 两种情况：
// 1. 当前i % pj == 0: pj为i的最小质因子，因此pj为pj * i的最小质因子
// 2. 当前i % pj == 1: pj一定小于i的所有质因子，因此pj一定为pj * i的最小质因子
// 每次均使用最小质因子筛，因此时间复杂度与n成正比
void get_primes2(int n) {
    for (int i = 2; i <= n; ++i) {
        if (!st[i]) primes[cnt++] = i;
        for (int j = 0; primes[j] <= n / i; ++j) {
            st[primes[j] * i] = 1;
            if (i % primes[j] == 0) break;
        }
    }
}

/////////////////////////////////////////////////////////

// 试除法求约数
/*
    时间复杂度分析：
    考虑1～n的所有数的约数
    1是所有数的约数，+n
    2，+n / 2
    ...
    所以1～n的约数个数：
    n / 1 + n / 2 + … + n / n = n(1 + 1 / 2 + … + n) = O(nlogn)
    单个数的约数：
    O(logn)
    因此排序时间复杂度：O(lognloglogn)
    对于单个数的总体时间复杂度：O(sqrt(n))
*/

// 计算所有约数（从小到大）
vector<int> get_factors(int x) {
    vector<int> ans;
    for (int i = 1; i <= x / i; ++i) {
        if (x % i == 0) {
            ans.push_back(i); // 较小的因数
            if (i != x / i) ans.push_back(x / i); // 较大的因数（相等只算一个）
        }
    }
    sort(ans.begin(), ans.end());
    return ans;
}

/////////////////////////////////////////////////////////

/*
    算数基本定理：
    N = p1 ^ a1 * ... * pn ^ an
    N的约数同样为p1～pn的组合，而ai的取值范围为[0, an]
    因此约数个数为(a1 + 1) * ... * (an + 1)

    约数之和即为所有质因数的组合的和
    (p1 ^ 0 + ... + p1 ^ a1) * ... * (pn ^ 0 + ... + pn ^ an)
*/

unordered_map<int, int> primes; // 质因数出现的次数
void get_prime_factors(int x) {
    for (int i = 2; i <= x / i; ++i) {
        while (x % i == 0) {
            x /= i;
            ++primes[i];
        }
    }
    if (x > 1) ++primes[x];
}

// 某个数的约数个数
LL get_factor_cnt() {
    LL res = 1;
    for (auto &p: primes) res = res * (p.second + 1) % MOD;
    return res;
}

// 某个数的约数之和
LL get_factor_sum() {
    LL res = 1;
    for (auto &[p, cnt]: primes) {
        LL t = 1;
        while (cnt--) t = (t * p + 1) % MOD;
        res = res * t % MOD;
    }
    return res;
}

// 1~n的约数之和（不是用算术基本定理，递推求解）
LL get_all_factor_sum(int n) {
    LL res = 0;
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j <= n; j += i) {
            ans += i;
        }
    }
    return ans;
}

/////////////////////////////////////////////////////////

// 欧几里得算法
// 求解最大公约数
int gcd(int a, int b) {
    return b ? gcd(b, a % b) : a;
}

// 扩展欧几里得算法
/*
    求解方程ax + by = gcd(a, b)
    使用gcd进行推广：
    (a, b) = (b, a % b);
    1. exgcd(a, 0, x, y) -> x = 1, y = 0
    2. exgcd(a, b, x, y) = exgcd(b, a % b, y, x) = d
       可以得到by + (a % b) * x
             = by + (a - (a / b) * b) * x
             = by + ax - (a / b) * bx
             = ax + b(y - (a / b) * x) = d
      -> x = x, y = y - (a / b) * x
    
    若计算出一组解为x_0, y_0，可以导出方程通解
    对于ax + by = gcd(a, b) = d
    令通解为x = (x_0 - k(b / d)), y = (y_0 + k(a / d))
    (代入可验证其满足ax_0 + by_0 = d)
*/
int exgcd(int a, int b, int &x, int &y) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    int d = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}

// 线性同余方程
/*
    求解：ax = b(mod m)
    可以将方程转化为：ax + my = b(mod m)
    然后求解d = exgcd(a, m, x, y)
    然后按照倍数关系计算x = b / d * x % m
    如果d不能整除b，则不存在解
*/
int d = exgcd(a, m, x, y);
if (b % d) x = (LL)b / d * x % m;

/////////////////////////////////////////////////////////

/*
    中国剩余定理（crt）
    方程组x = a_i(mod m_i)，其中1 <= i <= n，并且m_i两两互质
    构造解：
        令M = m_1 * ... * m_n, M_i = M / m_i，M_i_r为M_i的逆元(mod m_i)
        则x = a_1 * M_1 * M_1_r  + ... + a_n * M_n * M_n_r
    证明：
        将x代入第一个方程x = a_1(mod m_1)中，由于M_2, ..., M_n中均含有m_1，可以得到下式
        a_1 * M_1 * M_1_r = a_1(mod m_1)
        mod m_1时M_1 * M_1_r = 1
        因此得到a_1 = a_1(mod m_1)，对于第一个方程成立
        对于其余方程同理
    （由于m_i为素数，因此在计算M_i_r时可以采用pow（费马小定理要求是素数）或exgcd）
*/

int n;
LL a[N], m[N];

LL crt() {
    LL x = 0;
    LL M = 1;
    for (int i = 0; i < n; ++i) M *= m[i];
    for (int i = 0; i < n; ++i) {
        LL Mi = M / m[i];
        LL Mi_r, y;
        exgcd(Mi, m[i], Mi_r, y);
        // Mi_r = pow(Mi, m[i] - 2, m[i]);
        x = (x + a[i] * Mi * Mi_r % M) % M;
    }
    return (x % M + M) % M;
}

/*
    扩展中国剩余定理（excrt）：
    m_i并非两两互质
    构造递推式：
        对于第一个方程和第二个方程，进行变换
        x = a_1(mod m_1) -> x mod m_1 = a_1 -> x = k_1 * m_1 + a_1
        x = a_2(mod m_2) -> x mod m_2 = a_2 -> x = k_2 * m_2 + a_2
        于是有k_1 * m_1 + a_1 = k_2 * m_2 + a_2
        k_1 * m_1 - k_2 * m_2 = a_2 - a_1
        该方程有解等价于(m_1, m_2) | a_2 - a_1
        使用exgcd求解方程得到解k_1, k_2
        可以转化为通解形式k_1 + k * m_2 / d, k_2 + k * m_1 / d
        将第一个解代入第一个方程x = k_1 * m_1 + a_1中
        得到x = k * m_1 * m_2 / d + k_1 * m_1 + a_1
        对比观察可得
        a_1 = k_1 * m_1 + a_1
        m_1 = m_1 * m_2 / d
        最终a_1就是x的解（x = m + a(mod m) = a(mod m), 将更新后的a_1, m_1记为a, m）
*/

LL excrt() {
    LL m1 = m[0], a1 = a[0];
    for (int i = 1; i < n; ++i) {
        LL m2 = m[i], a2 = a[i];
        LL k1, k2;
        LL d = exgcd(m1, m2, k1, k2);
        if ((a2 - a1) % d) {
            a1 = -1; // -1表示无解
            break;
        }
        k1 *= (a2 - a1) / d;
        k1 = (k1 % (m2 / d) + (m2 / d)) % (m2 / d);

        a1 = k1 * m1 + a1;
        m1 = m1 / d * m2;
    }
    if (a1 != -1) a1 = (a1 % m1 + m1) % m1;
    return a1;
}

/////////////////////////////////////////////////////////

// 欧拉函数
/*
    [1, x]中与x互质的数的个数
    使用容斥原理计算：
    1. 总个数为x，首先减去包含因子pi的个数：x - x / p1 - ... - x / pn
    2. 加上包含因子pi * pj的个数：x - x / p1 - ... - x / pn + x / (p1 * p1) + ... + x / (pn * pn)
    ...
    合并得到phi(x) = x * (1 - 1 / p1) * ... * (1 - 1 / pn)
*/

int phi(int x) {
    int res = x;
    for (int i = 2; i <= x / i; ++i) {
        if (x % i == 0) {
            res = res / i * (i - 1);
            while (x % i == 0) x /= i;
        }
    }
    if (x > 1) res = res / x * (x - 1);
    return res;
}

/////////////////////////////////////////////////////////

// [1, n]的欧拉函数
/*
    1. 当x为质数时，phi(x) = x - 1
    2. 当i % pj == 0时（pj为i的最小质因子）
    phi(pj * i) = pj * i * ... * (1 - 1 / pj)
    phi(i) = i * ... * (1 - 1 / pj)
    因此phi(pj * i) = pj * phi(i)
    3. 当i % pj != 0时（pj小于i的所有质因子，pj为pj * i的最小质因子）
    phi(pj * i) = pj * i * ... * (1 - 1 / pj)
    phi(i) = i * ...
    因此phi(pj * i) = pj * phi(i) * (1 - 1 / pj) = (pj - 1) * phi(i)
*/

int primes[N], cnt;
int st[N];
int euler[N];

void get_euler(int n) {
    euler[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (!st[i]) {
            primes[cnt++] = i;
            euler[i] = i - 1;
        }
        for (int j = 0; primes[j] <= n / i; ++j) {
            int t = primes[j] * i;
            st[t] = 1;
            euler[t] = (primes[j] - 1) * euler[i];
            if (i % primes[j] == 0) {
                euler[t] = primes[j] * euler[i];
                break;
            }
        }
    }
}

/////////////////////////////////////////////////////////

// 高斯消元法
// 时间复杂度：O(n^3)

int n;
double a[N][N];

int gauss_elimination() {
    int c, r;
    for (c = 0, r = 0; c < n; ++c) {
        // 找到当前列的值里绝对值最大的一行
        int t = r;
        for (int i = r; i < n; ++i) {
            if (fabs(a[i][c]) > fabs(a[t][c])) t = i;
        }
        if (fabs(a[t][c]) < eps) continue;
        for (int i = c; i <= n; ++i) swap(a[t][i], a[r][i]); // 将该行换到第r行
        for (int i = n; i >= c; --i) a[r][i] /= a[r][c]; // 使该行第一个元素为1
        for (int i = r + 1; i < n; ++i) { // 用当前行把所在的第一列全消为0
            if (fabs(a[i][c]) > eps) {
                for (int j = n; j >= c; --j) {
                    a[i][j] -= a[r][j] * a[i][c];
                }
            }
        }
        ++r;
    }
    if (r < n) {
        for (int i = r; i < n; ++i) {
            if (fabs(a[i][n]) > eps) return 2; // 无解
        }
        return 1; // 无穷多解
    }
    // 把对角线上方元素消为0（实际效果是直接将结果作用到最后一列）
    // e.g.
    // 1 2 3 4          1 2 3 4          1 2 0 1          1 0 0 -1
    // 0 1 2 3-(3-1*2)->0 1 0 1-(4-1*3)->0 1 0 1-(1-1*2)->0 1 0 1
    // 0 0 1 1          0 0 1 1          0 0 1 1          0 0 1 1
    for (int i = n - 2; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            a[i][n] -= a[i][j] * a[j][n];
        }
    }
    return 0; // 唯一解
}

/////////////////////////////////////////////////////////

// 组合数（暴力递推）
// 时间复杂度：O(n^2)
// c[i][j] = c[i - 1][j - 1] + c[i - 1][j]

void init() {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            if (!j) c[i][j] = 1;
            else c[i][j] = (c[i - 1][j - 1] + c[i - 1][j]) % MOD;
        }
    }
}

// 组合数（预处理阶乘逆元）大量常用
// 时间复杂度：O(nlogn)
// c[a][b] = a! / (b! * (a - b)!)

void init() {
    fact[0] = infact[0] = 1;
    for (int i = 1; i < N; ++i) {
        fact[i] = (LL)fact[i - 1] * i % MOD;
        infact[i] = (LL)infact[i - 1] * pow(i, MOD - 2, MOD) % MOD;
    }
}

(LL)fact[a] * infact[b] % MOD * infact[a - b] % MOD; // c[a][b]

// 组合数（卢卡斯定理递推）
// 单次时间复杂度：O(nlogn)
/*
    卢卡斯定理：
    c[a][b] = c[a mod p][b mod p] * c[a / p][b / p]，其中p为质数
    证明：
    将a表示为k进制
    a = a_k * p ^ k + ... + a_0 * p ^ 0
    则(1 + x) ^ a = ((1 + x) ^ (p ^ k)) ^ a_k * ... * ((1 + x) ^ (p ^ 0)) ^ a_0
    而(1 + x) ^ p = 1 + c[p][1] * x ^ 1 + ... + c[p][p] * x ^ p(mod p)
                  = 1 + x ^ p(mod p)
    因此(1 + x) ^ a = (1 + x ^ (p ^ k)) ^ a_k * ... * (1 + x ^ (p ^ 0)) ^ a_0(mod p)
    注意到b = b_k * p ^ k + ... + b_0 * p ^ 0
    因此x ^ b的系数为c[a][b] = c[a_k][b_k] * ... * c[a_0][b_0]
    而c[a][b] = c[a mod p][b mod p] * c[a / p][b / p]
              = c[a mod p][b mod p] * c[(a / p) mod p][(b / p) mod p] * c[(a / p) / p][(b / p) / p]
              ...
              = c[a_k][b_k] * ... * c[a_0][b_0]
*/

int C(LL a, LL b, int p) {
    if (b > a) return 0;
    int res = 1;
    for (int i = 1, j = a; i <= b; ++i, --j) {
        res = (LL)res * j % p;
        res = (LL)res * pow(i, p - 2, p) % p;
    }
    return res;
}

// c[a][b]
int lucas(LL a, LL b, int p) {
    if (a < p && b < p) return C(a, b, p);
    return (LL)C(a % p, b % p, p) * lucas(a / p, b / p, p) % p;
}

// 组合数（线性求逆元）单个常用
// c[n][i] = (n - i + 1) / i * c[n][i - 1]
// 时间复杂度O(n)
/*
    线性求逆元i ^ (-1)
    推导：
    1. 当i = 1时，i ^ (-1) = 1(mod p)
    2. 当i > 1时
       令k = p / i, j = p mod i, 得p = k * i + j
       k * i + j = 0(mod p)
       两边同时乘上i ^ (-1) * j ^ (-1)
       k * j ^ (-1) + i ^ (-1) = 0(mod p)
       i ^ (-1) = -k * j ^ (-1)(mod p)
       即i ^ (-1) = -(p / i) * (p mod i) ^ (-1)(mod p)
       由于p mod i < i, 因此可以使用(p mod i) ^ (-1)推导i ^ (-1)
       为了防止负数，-(p / i)加上p
       i ^ (-1) = (p - p / i) * (p mod i) ^ (-1)(mod p)
*/

int inv[N];

inv[1] = 1;
for (int i = 2; i <= n; ++i) inv[i] = (LL)(MOD - MOD / i) * inv[MOD % i] % MOD;

// 线性递推求组合数：c[n][i] = (n - i + 1) / i * c[n][i - 1]
for (int i = 1; i <= k; ++i) c[n][i] = (LL)c[n][i - 1] * (n - i + 1) % MOD * inv[i] % MOD;

/////////////////////////////////////////////////////////

// 卡特兰数
/*
    Catalan(n) = 1 / (n + 1) * c[2n][n]
    证明:
    使用画图即可（从(0, 0)到(n, n)的路径）
    将第一次经过非法边界以后到终点的路线以非法边界为轴对称
    则合法的结果就是所有结果减去非法结果
    c[2n][n] - c[2n][n + 1] = 1 / (n + 1) * c[2n][n]
*/

(LL)c[2 * n][n] * inv[n + 1] % MOD

/////////////////////////////////////////////////////////

/*
    容斥原理：
    |S_1 + ... + S_n| = ∑|S_i| - ∑|S_i + S_j| + ∑|S_i + S_j + S_k| ...
    2 ^ n - 1 = c[n][1] + ... + c[n][n]
*/

/////////////////////////////////////////////////////////

// 简单博弈
/*
    Nim游戏：（以下^为异或）
    1. 拿完：0 ^ ... ^ 0 = 0（必败）
    2. a_1 ^ ... ^ a_n = x != 0（必胜）
       假设x的最高位1在第k位
       a_1, ... , a_n中必有一个a_i的第k位为1
       则a_i ^ x < a_i
       可以从第a_i中拿走a_i - a_i ^ x个
       则a_1 ^ ... ^ (a_i ^ x) ^ ... ^ a_n = x ^ x = 0
       因此可以由当前状态转为必败状态，当前为必胜态
    3. a_1 ^ ... ^ a_n = x = 0（必败）
       假设存在一种方法，使得a_1 ^ ... ^ a_i' ^ ... ^ a_n = 0 (a_i != a_i')
       则a_1 ^ ... ^ a_i ^ ... ^ a_n ^ a_1 ^ ... ^ a_i ^ ... ^ a_n = 0
       a_i ^ a_i' = 0
       a_i = a_i', 与a_i != a_i'矛盾
       因此不存在一种操作将当前状态变为必败态，当前为必败态
*/

/*
    台阶Nim游戏：把台阶上的石头往下拿
    因为第一阶拿到地面要拿一次，第二阶拿两次，第三阶拿三次…
    所以可以看成第二阶有两堆石子，第三阶有三堆....
    因为偶数阶石子为偶数堆，所以异或为0，奇数阶异或后就是原本石子数目
    所以可以把原本所有奇数阶的石子进行异或，得到的就是答案
*/

/*
    SG定理：
    定义mex函数的值为不属于集合S中的最小非负整数
    SG(x) = mex{SG(y_1), ..., SG(y_k)}
    y_i为x可以到达的状态，终点的SG为0
    如果有状态SG(x_1), ... , SG(x_n)
    则SG(x_1) ^ ... ^ SG(x_n) != 0必胜
    证明同从几个堆中选定一个然后从中取石子类似
*/

int k;
int s[N]; // 可以取的k个种类
int f[M];

// 记忆化搜索
int sg(int x) {
    if (f[x] != -1) return f[x];
    unordered_set<int> st;

    // 每次只能取s[i]个
    for (int i = 0; i < k; ++i) {
        if (x >= s[i]) st.insert(sg(x - s[i]));
    }
    // 每一堆可以替换成比当前小的2堆
    // for (int i = 0; i < x; ++i) {
    //     for (int j = 0; j <= i; ++j) {
    //         st.insert(sg(i) ^ sg(j));
    //     }
    // }

    for (int i = 0;; ++i) {
        if (!st.count(i)) {
            return f[x] = i;
        }
    }
}

/////////////////////////////////////////////////////////

// 背包问题

/*
    01背包：
    状态表示：
        f[i][j]: 
            集合：只考虑前i个，体积不大于j
            属性：max
    状态计算：
        集合划分：
            不含i：f[i - 1][j]
            含i：f[i - 1][j - v[i]] + w[i]

    由于f[i][j] = max(f[i - 1][j], f[i - 1][j - v[i]] + w[i])
    因此考虑使用滚动数组优化
    f[j] = max(f[j], f[j - v[i]] + w[i])
    但是此时j - v[i] < j，对于当前i在更新f[j]时，如果f[j - v[i]]先更新
    则使用的实际上是f[i][j - v[i]]
    因此j应该反向遍历
*/

int n, m;
int v[N], w[N];
int f[N];

for (int i = 1; i <= n; ++i) {
    for (int j = m; j >= v[i]; --j) {
        f[j] = max(f[j], f[j - v[i]] + w[i]);
    }
}

/*
    01背包的初始化
    对于体积j：
        至少为j：j<=0初始化为0，其余为不合法
        恰好为j：j=0初始化为0，其余为不合法
        至多为j：全部初始化为0
*/

/*
    完全背包：
    状态表示：
        f[i][j]: 
            集合：只考虑前i个，提及不大于j
            属性：max
    状态计算：
        集合划分：
            不含i：f[i - 1][j]
            含k个i：f[i - 1][j - k * v[i]] + k * w[i]

    f[i][j] = max(f[i - 1][j], f[i - 1][j - v[i]] + w[i], f[i - 1][j - 2 * v[i]] + 2 * w[i], ...)
    f[i][j - v[i]] = max(f[i - 1][j - v[i]], f[i - 1][j - 2 * v[i]] + w[i], ...)
    观察对比可知f[i][j] = max(f[i - 1][j], f[i][j - v[i]] + w[i])
    因此考虑使用滚动数组优化
    f[j] = max(f[j], f[j - v[i]] + w[i])
    此时j - v[i] < j，f[j - v[i]]实际上是当前i的f[i][j - v[i]]
    因此j应该正向遍历
*/

int n, m;
int v[N], w[N];
int f[N];

for (int i = 1; i <= n; ++i) {
    for (int j = v[i]; j <= m; ++j) {
        f[j] = max(f[j], f[j - v[i]] + w[i]);
    }
}

/*
    多重背包（暴力）
    时间复杂度：O(NVS)
    f[i][j] = max(f[i - 1][j - k * v[i]] + k * w[i]), k = 0, ..., s[i]
*/

for (int i = 1; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
        for (int k = 0; k <= s[i] && j - k * v[i] >= 0; ++k) {
            f[i][j] = max(f[i][j], f[i - 1][j - k * v[i]] + k * w[i]);
        }
    }
}

/*
    多重背包（二进制优化）
    时间复杂度：O(NVlogS)
    f[i][j] = max(f[i - 1][j], f[i - 1][j - v[i]] + w[i], ... , f[i - 1][j - s[i]v[i]] + s[i]w[i])
    f[i][j - v[i]] = max(f[i - 1][j - v[i]], ... , f[i - 1][j - (s[i] + 1) * v[i]] + s[i] * w[i])
    f[i][j - v[i]]多出一项f[i - 1][j - (s[i] + 1) * v[i]] + s[i] * w[i]
    因此并不能采用代入的优化方法
    从f[i][j]到f[i][j - v[i]]，相当于需要计算滑动窗口最大值
    将所有s[i]分解成可以凑成[1, s[i]]的二的倍数
    e.g. s[i] = 9, 使用1, 2, 4, 2可以凑出[1, 9]
    1, 2, 4 -> 1, 2, 3, 4, 5, 6, 7
    多一个2 -> 3, 4, 5, 6, 7, 8, 9
    => 1, 2, 3, 4, 5, 6, 7, 8, 9
    一般地，使用1, 2, 4, ..., x可以凑出来[1, 2 * x - 1]
    加上一个剩下部分可以凑出所有需要的数字
    然后重新计算v[i]和w[i]之后放入到v和w数组中
    问题转化为01背包
    时间复杂度从O(NVS)降到O(NVlogS)
*/

// N取N * logS
int n, m;
int v[N], w[N];
int f[M];

// 二进制预处理
int cnt = 1;
for (int i = 1; i <= n; ++i) {
    int vv, ww, ss;
    cin >> vv >> ww >> ss;
    int k = 1;
    while (k <= ss) {
        v[cnt] = vv * k;
        w[cnt] = ww * k;
        ++cnt;
        ss -= k;
        k <<= 1;
    }
    if (ss) {
        v[cnt] = vv * ss;
        w[cnt] = ww * ss;
        ++cnt;
    }
}
// 转化成01背包
for (int i = 1; i <= cnt; ++i) {
    for (int j = m; j >= v[i]; --j) {
        f[j] = max(f[j], f[j - v[i]] + w[i]);
    }
}

/*
    多重背包（单调队列优化）
    时间复杂度：O(NV)
    对于权重w的处理：
    当前队列中元素的范围是[r, j]，当前元素是j
    下标减小减少权重，下标增大增加权重
    f[i][j] = max(f[i - 1][j], ... , f[i - 1][j - s * v] + s * w)
    假设当前下标为j，队头下标为j - s * v
    未修正的比较：f[j] v.s. f[j - s * v]
    修正后的比较：f[j] - (j - r) / v * w v.s. f[j - s * v] - (j - s * v - r) / v * w（均修正到r进行比较）
    即（两边同时加上(j - r) / v * w）: f[j] v.s. f[i - 1][j - s * v] + s * w
    与状态方程中表达的一致
*/

#include <iostream>
#include <cstring>

using namespace std;

const int N = 100010;

int n, m;
int f[N], g[N], q[N];

for (int i = 0; i < n; ++i) {
    int v, w, s;
    cin >> v >> w >> s;
    memcpy(g, f, sizeof f); // 将f复制到g
    for (int r = 0; r < v; ++r) { // 枚举可能的余数
        int hh = 0, tt = -1;
        for (int j = r; j <= m; j += v) { // 枚举当前体积
            while (hh <= tt && q[hh] < j - s * v) ++hh; // 滑出窗口
            while (hh <= tt && g[q[tt]] - (q[tt] - r) / v * w <= g[j] - (j - r) / v * w) --tt; // 队尾比当前小，队尾出队
            q[++tt] = j;
            f[j] = g[q[hh]] + (j - q[hh]) / v * w; // 窗口最大值即为队头
        }
    }
}

/*
    分组背包
    每个组选一个，和01背包类似，多了一个当前组内选择哪一个的维度
*/

int n, m;
int s[N], v[N][N], w[N][N];
int f[N][N];

for (int i = 1; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
        f[i][j] = f[i - 1][j];
        for (int k = 1; k <= s[i]; ++k) {
            if (j >= v[i][k])
                f[i][j] = max(f[i][j], f[i - 1][j - v[i][k]] + w[i][k]);
        }
    }
}

/*
    分组背包（优化空间）
*/

int n, m;
int s[N], v[N][N], w[N][N];
int f[N];

for (int i = 1; i <= n; ++i) {
    for (int j = m; j >= 0; --j) {
        for (int k = 0; k < s[i]; ++k) {
            if (j >= v[i][k]) f[j] = max(f[j], f[j - v[i][k]] + w[i][k]);
        }
    }
}

/*
    二维费用背包
    需要同时考虑体积、质量
    类似01背包
*/

int n, V, M;
int f[N][N];

for (int i = 0; i < n; ++i) {
    int v, m, w;
    cin >> v >> m >> w;
    for (int j = V; j >= v; --j) {
        for (int k = M; k >= m; --k) {
            f[j][k] = max(f[j][k], f[j - v][k - m] + w);
        }
    }
}

/*
    01背包求具体方案
    使用非空间优化的版本方便回退
*/

int g[N]; // 逆序答案

// 回退计算答案
int j = m;
for (int i = n; i >= 1; --i) {
    for (int k = 1; k <= m; ++k) {
        if (j >= v[i][k] && f[i][j] == f[i - 1][j - v[i][k]] + w[i][k]) {
            j -= v[i][k];
            g[i] = k;
            break;
        }
    }
}
// 答案即为g[1:n]










/////////////////////////////////////////////////////////

// 最长上升子序列（LIS）
/*
    d[i]: 长度为i的最小结尾元素
    如果a[i] > d[len]，则长度可以加1
    否则找到第一个大于等于a[i]的元素进行更新（二分）
    时间复杂度：O(nlogn)
*/

int a[N], d[N];

d[++len] = a[0];
for (int i = 0; i < n; ++i) {
    if (a[i] > d[len]) d[++len] = a[i];
    else {
        int l = 1, r = len;
        while (l < r) {
            int mid = (l + r) >> 1;
            if (d[mid] >= a[i]) r = mid;
            else l = mid + 1;
        }
        d[l] = a[i];
    }
}

/////////////////////////////////////////////////////////

// 最长公共子序列（LCS）
/*
    f[i][j]: 当前a[1:i]和b[1:j]的最长公共子序列
    子集划分为：00 01 10 11
    分别表示当前位置i, j选或不选
    00: 使用f[i - 1][j - 1]表示
    01, 10: 使用f[i - 1][j], f[i][j - 1]来转移
    11: 使用f[i - 1][j - 1]来转移

    使用以上方法得到的4个子集可能有交集
    例如01一定包含j，而f[i - 1][j]不一定包含j
    因此f[i - 1][j]一定包含01
    由于计算的是max，因此并不会影响最终结果

    f[i - 1][j - 1]已经包含在f[i - 1][j], f[i][j - 1]中
    因此只需要计算f[i - 1][j], f[i][j - 1], f[i - 1][j - 1] + 1
*/

int n, m;
char a[N], b[N];
int f[N][N];

for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
        f[i][j] = max(f[i - 1][j], f[i][j - 1]);
        if (a[i] == b[j]) f[i][j] = max(f[i][j], f[i - 1][j - 1] + 1);
    }
}

/////////////////////////////////////////////////////////

// 编辑距离
// 增删改A得到B的最少次数

int n, m;
char a[N], b[N];
int f[N][N];

for (int i = 0; i <= n; ++i) f[i][0] = i;
for (int j = 0; j <= m; ++j) f[0][j] = j;
for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
        if (a[i] == b[j]) f[i][j] = f[i - 1][j - 1];
        else f[i][j] = min({f[i - 1][j - 1], f[i - 1][j], f[i][j - 1]}) + 1;
    }
}

/////////////////////////////////////////////////////////

// 区间DP
/*
    将当前区间[i, j]划分为[i, k], [k + 1, r]:
        先枚举区间长度len
        再枚举区间起点i
        再枚举区间分割点k
    时间复杂度：O(n ^ 3)
*/

for (int len = 2; len <= n; ++len) {
    for (int i = 1; i + len - 1 <= n; ++i) {
        int l = i, r = i + len - 1;
        for (int k = l; k < r; ++k) {
            // f[l][r] = ...
        }
    }
}

/////////////////////////////////////////////////////////

// 状态压缩DP

/*
    1.棋盘型（连通型）
        先枚举转移
        后枚举状态
    2.集合型
        先枚举状态
        后枚举转移
*/

// 在一个n*m矩形中放1*2的砖块的可能方案数（棋盘型，连通型）
/*
    该问题只需要确定横向砖块的摆放，即可得到纵向砖块的唯一摆法
    f[i][j]: i为当前列号, j为当前列的状态（二进制）
    当前位置为1：1*2砖块的右半部分
    合法摆法的相邻两列状态j，k需要满足两个条件：
    (j & k) == 0
    (j | k)不存在连续奇数个0（连续偶数个0可以竖放）
    合法的最终状态相当于第m + 1列什么都不放
*/

// 最短哈密顿路径（集合型）
/*
    f[i][j]: 所有0到j路径中包含状态i（二进制）的最短路
    f[i][j] = min(f[i][j], f[i - {j}][k] + w[k][j])
    其中k为直连j的上一个点
*/

/////////////////////////////////////////////////////////

// 树形DP

void dfs(int u) {
    // f[u][0] = ...
    // f[u][1] = ...
    for (int i = h[u]; ~i; i = ne[i]) {
        int v = e[i];
        dfs(v);
        // f[u][0] = ...
        // f[u][1] = ...
    }
}

// 树中最长路径
/*
    该长度一定是某个点往下延申的最大与次大长度之和
    对于某个点u，其所有子节点为v（不能包含u的父节点p）
    用dfs(v, u) + w[u->v]更新最大值与次大值
    用最大值与次大值之和更新答案
*/
int dfs(int u, int p) { // u：当前，p：父
    int d1 = 0, d2 = 0; // 最大，次大

    for (int i = h[u]; ~i; i = ne[i]) {
        int v = e[i];
        if (v == p) continue;
        int d = dfs(v, u) + w[i];
        if (d >= d1) d2 = d1, d1 = d;
        else if (d > d2) d2 = d;
    }
    ans = max(ans, d1 + d2);
    return d1;
}

/////////////////////////////////////////////////////////

// 记忆化搜索

int f[N];

memset(f, -1, sizeof f);

int dfs(int u) {
    int &res = f[u];
    if (res != -1) return res;
    // res = ...
    for (int i = h[u]; ~i; i = ne[i]) {
        int v = e[i];
        // res = ...
    }
    f[u] = res;
    return res;
}

/////////////////////////////////////////////////////////

/*
    迭代式数位DP：
    f[i][j][*]: 一共i位，最高位为j，状态为*的函数（可以是一个变量、或者一个包含多个变量的结构体）
    以下为基本框架（具体情况具体分析）
*/

void init() {
    for (int j = 0; j <= 9; ++j) f[1][j][*] = *;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= 9; ++j) {
            for (int k = 0; k <= 9; ++k) {
                f[i][j][*] <= f[i - 1][k][*]; // 转移计算
            }
        }
    }
}

// 计算时候采用前缀和的思想，dp[l:r] = dp(r) - dp(l - 1)
int dp(int n) {
    if (!n) return *;
    vector<int> nums;
    while (n) nums.push_back(n % 10), n /= 10;
    int ans = 0;
    int last = *; // 上一个数字（例如：上一个数、之前所有数的和）
    for (int i = nums.size() - 1; i >= 0; --i) {
        int x = nums[i];
        for (int j = 0; j < x; ++j) { // 当前位为0~x-1，后面可以随便填，累加计算合法的数量
            ans += f[i + 1][j];
        }
        if (x不满足某种条件) break;
        last <= x; // 更新last
        if (!i && last满足某种条件) ++ans; // 整体满足某种条件，答案加一
    }
    return ans;
}

/////////////////////////////////////////////////////////

/*
    哈夫曼树证明：（不严谨）
    最小的两个点必然是深度最深并且可以是兄弟
    如果不是，则可以将深度小的点进行替换，得到更优的答案
    因此得到最优解的第一步，必然可以是合并最小的两个点
    将最小2点合并后，规模从n变为n - 1，继续通过合并当前最小的2个点
    最终得到最优解
*/
// 使用优先级队列每次拿出2个合并

/////////////////////////////////////////////////////////

/*
    快速傅里叶变换（FFT）
    对于C(x) = A(x)B(x)的计算
    1. 加倍：将A，B变为2n次（高次补0）
    2. 求值：2n阶FFT计算A与B的点值表达
    3. 逐点相乘：C = AB
    4. 插值：C的系数c = DFT^{-1}(C)

    DFT
    n次单位复根：满足w^n = 1的w，恰好n个(e^{2*pi*i*k/n}, k=0,...,n-1)
    记n次主单位根为w_n = e^{2*pi*i/n}, 则其他所有n次单位复根都是w_n的幂次

    消去引理：w_{dn}^{dk} = w_{n}^{k}, n>=0, k>=0, d>0
    推论：对于任意偶数n>0，有w_{n}^{n/2} = w_{2} = -1

    折半引理：对于偶数n>0，n个n次单位复根平方的集合等价于n/2个n/2次单位复根的集合
    证明：(w_{n}^{k+n/2})^2 = (w_{n}^{k})^2
    即w_{n/2}^{k+n/2} = w_{n/2}^{k}
    该集合每2个这样的元素是相等的，因此集合大小为n/2

    求和引理：n>=1, k不能被n整除
    sum_{j=0}^{n-1}((w_{n}^{k})^j) = 0

    A(x) = sum(a_j*x^j)
    y_k = A(w_{n}^{k})
    y = DFT_{n}(a), k=0,...,n-1

    令A0(x) = a0 + a2*x + ... + a_{n-2}x^{n/2-1}
      A1(x) = a1 + a3*x + ... + a_{n-1}x^{n/2-1}
    有A(x) = A0(x^2) + x * A1(x^2)
    分别计算A0, A1在(w_{n}^{0})^2, ... , (w_{n}^{n-1})^2处的值，综合得到A

    以上过程可以使用FFT来实现
    FFT(a):
    n = a.size()
    if n == 1
        return a
    w_n = e^{2*pi*i/n}
    w = 1
    a0 = (a_0, a_2, ..., a_{n-2})
    a1 = (a_1, a_3, ..., a_{n-1})
    y0 = FFT(a0)
    y1 = FFT(a1)
    for k = 0 to n/2 - 1
        t = w * y1_{k}
        u = y0_{k}
        y_{k} = u + t
        y_{k+n/2} = u - t
        w = w * w_n
    return y

    接下来剩下最后一步c = DFT^{-1}(C)
    将DFT写为：y = V_{n} * a
    则a = V_{n}^{-1} * y
    可以证明V_{n}^{-1}的(j, k)处元素为(w_{n}^{-k*j}) / n
    因此将FFT修改可以得到逆FFT
    1. a与y互换
    2. 使用w_{n}^{-1}替换w_{n}
    3. 结果的每个元素除以n

    多项式乘法可以表示为卷积定理的形式
    卷积定理：a x b = DFT^{-1}(DFT(a) * DFT(b)), 式中a, b均为2n阶
    
    注意到递归FFT过程中的计算顺序，e.g. n=8
    (a0, a1, a2, a3, a4, a5, a6, a7)
    (a0, a2, a4, a6)(a1, a3, a5, a7)
    (a0, a4)(a2, a6)(a1, a5)(a3, a7)
    (a0)(a4)(a2)(a6)(a1)(a5)(a3)(a7)
    可以直接使用(a0)(a4)(a2)(a6)(a1)(a5)(a3)(a7)的顺序迭代计算FFT

    此处需要一个change算法将y的顺序调整到最终状态
    从小到大递推实现
    设len=2^k，其中k为二进制数的长度，设R(x)为二进制数x反转后的数
    显然，R(0)=0
    从小到大计算，因此在计算R(x)时，R(x/2)的值已知
    因此将x右移一位，然后翻转，再右移一位，最后在结果最高位加上x最低位的值
    即R(x) = R(x/2)/2 + (x%2) * (len/2)

    最终得到迭代FFT
    FFT(a):
    change(a)
    n = n = a.size()
    for s = 1 to log{n}
        m = 2^s
        w_m = e^{2*pi*i/m}
        for k = 0 to (n - 1) by m
            w = 1
            for j = 0 to m/2 - 1
                t = w * a[k + j + m/2]
                u = a[k + j]
                a[k + j] = u + t
                a[k + j + m/2] = u - t
                w = w * w_m
    return a
*/

// 使用时候要注意N的范围，扩展到大于等于指定长度的2的次幂有可能会超过N
// 通常开到3倍

#include <cmath>

using namespace std;

const int N = 100010;

const double PI = acos(-1.0);

struct Complex {
    double x, y;
    Complex(double _x = 0.0, double _y = 0.0) : x(_x), y(_y) {}
    Complex operator-(const Complex &rhs) const {
        return Complex(x - rhs.x, y - rhs.y);
    }
    Complex operator+(const Complex &rhs) const {
        return Complex(x + rhs.x, y + rhs.y);
    }
    Complex operator*(const Complex &rhs) const {
        return Complex(x * rhs.x - y * rhs.y, x * rhs.y + y * rhs.x);
    }
};

int rev[N];
Complex x1[N * 2], x2[N * 2];
char str1[N], str2[N];
int sum[N * 2];

void change(Complex y[], int n) {
    for (int i = 0; i < n; ++i) rev[i] = i;
    for (int i = 0; i < n; ++i) {
        rev[i] = rev[i >> 1] >> 1;
        if (i & 1) rev[i] |= n >> 1;
    }
    for (int i = 0; i < n; ++i) {
        if (i < rev[i]) swap(y[i], y[rev[i]]);
    }
}

// fft: r=1, ifft: r=-1
void fft(Complex y[], int n, int r) {
    change(y, n);
    for (int m = 2; m <= n; m <<= 1) {
        Complex w_m(cos(2 * PI / m), sin(r * 2 * PI / m));
        for (int k = 0; k < n; k += m) {
            Complex w(1, 0);
            for (int j = 0; j < m / 2; ++j) {
                Complex t = w * y[k + j + m / 2];
                Complex u = y[k + j];
                y[k + j] = u + t;
                y[k + j + m / 2] = u - t;
                w = w * w_m;
            }
        }
    }
    if (r == -1) {
        for (int i = 0; i < n; ++i) y[i].x /= n;
    }
}

// 以下计算2个大整数str1, str2乘积
int l1 = strlen(str1), l2 = strlen(str2);
int l = 1;
while (l < l1 * 2 || l < l2 * 2) l <<= 1; // 长度扩展到不小于l1, l2的2的倍数
for (int i = 0; i < l1; ++i) x1[i] = Complex(str1[l1 - i - 1] - '0', 0);
for (int i = l1; i < l; ++i) x1[i] = Complex(0, 0);
for (int i = 0; i < l2; ++i) x2[i] = Complex(str2[l2 - i - 1] - '0', 0);
for (int i = l2; i < l; ++i) x2[i] = Complex(0, 0);
fft(x1, l, 1); // 求值
fft(x2, l, 1); // 插值
for (int i = 0; i < l; ++i) x1[i] = x1[i] * x2[i];
fft(x1, l, -1); // 计算
for (int i = 0; i < l; ++i) sum[i] = int(x1[i].x + 0.5); // 四舍五入
for (int i = 0; i < l; ++i) { // 进位
    sum[i + 1] += sum[i] / 10;
    sum[i] %= 10;
}
l = l1 + l2 - 1;
while (l && sum[l] == 0) --l; // 删除前导0
for (int i = l; i >= 0; --i) printf("%d", sum[i]);

/////////////////////////////////////////////////////////


// 矩阵旋转
/*
    (x, y)
    (y, n - x - 1) // 顺时针旋转90
    (n - y - 1, x) // 逆时针旋转90
    更高度数依此类推
*/
    
/////////////////////////////////////////////////////////