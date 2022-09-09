import scipy
from scipy.io import arff

k = int(input('Enter k: '))
input_file = input('Enter arff file: ')

data = arff.loadarff(input_file)[0]

class Graph:
    def __init__(self, n):
        self.n = n
        self.edge = []
        self.mst = None

    def insert_edge(self, u, v, w):
        self.edge.append([u, v, w])

    # union find
    def find(self, par, i):
        if par[i]==i:
            return i
        par[i] = self.find(par, par[i])
        return self.find(par, par[i])

    # union function for union-find structure
    def union(self, par, rank, u, v):
        ur = self.find(par, u)
        vr = self.find(par, v)

        if (rank[ur] < rank[vr]):
            par[ur] = vr
        elif (rank[ur] > rank[vr]):
            par[vr] = ur
        else:
            par[vr] = ur
            rank[ur] += 1

    # kruskal
    def MST(self):
        mst, par, rank = [], [], []
        for i in range(self.n):
            par.append(i)
            rank.append(0)
        i, e = 0, 0
        self.edge = sorted(self.edge, key = lambda edg: edg[2])
        while e < self.n - 1:
            u, v, w = self.edge[i]
            i = i+1
            x, y = self.find(par, u), self.find(par, v)
            if x!=y:
                e = e+1
                mst.append([u, v, w])
                self.union(par, rank, x, y)
        return mst

    # k-clustering
    def clustering(self, k):
        if self.mst == None:
            self.mst = self.MST()
        visited = [-1 for i in range(self.n)]
        adj = [[] for i in range(self.n)]
        # removing k-1 expensive edges
        for i in range(0, self.n - k):
            adj[self.mst[i][0]].append(self.mst[i][1])
            adj[self.mst[i][1]].append(self.mst[i][0])
        cur_cluster = 0
        # bfs to find connected components in the forest
        for i in range(self.n):
            if (visited[i] == -1):
                queue = [i]
                visited[i] = cur_cluster
                while len(queue) != 0:
                    u = queue.pop(0)
                    for v in adj[u]:
                        if visited[v] == -1:
                            visited[v] = cur_cluster
                            queue.append(v)
                cur_cluster = cur_cluster + 1
        clusters = [[] for i in range(cur_cluster)]
        for i in range(self.n):
            clusters[visited[i]].append(i)
        return clusters

# function to calculate d
def d(u, v):
    diff = 0
    for i in range(len(data[0])-1):
        diff = diff + (data[u][i] - data[v][i])*(data[u][i] - data[v][i])
    return diff**0.5

data = data[0:2000] 
g = Graph(len(data))
for i in range(g.n):
    for j in range(0, i):
        if (i!=j):
            g.insert_edge(i, j, d(i, j))

clusters = g.clustering(k)
print("START")
print(clusters)
print("END")

# calculating the purity
val = 0
for i in range(k):
    d = {}
    for u in clusters[i]:
        d[data[u][-1]] = d.get(data[u][-1], 0) + 1
    val = val + max(d.values())
print("Purity for " + str(k) + " clusters is: " + str(val / g.n))



purities = []
max_purity = 0
optimal_k = 0
for j in range(1, g.n + 1):
    val = 0
    clusters = g.clustering(j)
    for i in range(j):
        d = {}
        for u in clusters[i]:
            d[data[u][-1]] = d.get(data[u][-1], 0) + 1
        val = val + max(d.values())
        
    purity = val / g.n
    purity_updated = purity * (g.n - j) / g.n # penalizes very high purities
    if (purity_updated > max_purity):
        max_purity = purity_updated
        optimal_k = j
    purities.append(purity) 

print('max purity: ' + str(max_purity))
print('optimal k: ' + str(optimal_k))

running_val = 0 
optimal_clusters = g.clustering(optimal_k)
vals = {} 
for i in range(optimal_k): 
    d = {} 
    for u in optimal_clusters[i]: 
        d[data[u][-1]] = d.get(data[u][-1], 0) + 1
    val = max(d.values())
    most_prevalent_label = max(d, key=d.get)
    vals[i + 1] = (most_prevalent_label, val)
    running_val = running_val + val

purity = running_val / g.n
for cluster in range(1, optimal_k + 1):
    print("Most prevalent label for cluster " + str(cluster) + " is: " + str(vals[cluster][0]))

print("Purity is: " + str(purity))
print("Formula for calculating purity:")
numerator = ""
for cluster in range(1, optimal_k + 1):
    numerator = numerator + str(vals[cluster][1]) + "+"
print("1/" + str(g.n) + "*(" + numerator + ")")


import matplotlib.pyplot as plt
plt.plot(purities)
plt.xlabel("k")
plt.ylabel("Purity")
plt.savefig("plot.png")

res = ""
for i in range(g.n):
    if i%10 == 0:
        res = res + ", " + str(purities[i])

