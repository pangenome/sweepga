/// Union-Find (Disjoint Sets) data structure for transitive mapping merging
pub struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    /// Create a new UnionFind with n elements
    pub fn new(n: usize) -> Self {
        let parent = (0..n).collect();
        let rank = vec![0; n];
        UnionFind { parent, rank }
    }

    /// Find the root of element x with path compression
    pub fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    /// Union two sets containing x and y
    pub fn union(&mut self, x: usize, y: usize) {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x != root_y {
            // Union by rank
            if self.rank[root_x] < self.rank[root_y] {
                self.parent[root_x] = root_y;
            } else if self.rank[root_x] > self.rank[root_y] {
                self.parent[root_y] = root_x;
            } else {
                self.parent[root_y] = root_x;
                self.rank[root_x] += 1;
            }
        }
    }

    /// Check if two elements are in the same set
    pub fn connected(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }

    /// Get all sets as groups of indices
    pub fn get_sets(&mut self) -> Vec<Vec<usize>> {
        let n = self.parent.len();
        let mut root_to_group: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();

        for i in 0..n {
            let root = self.find(i);
            root_to_group.entry(root).or_insert_with(Vec::new).push(i);
        }

        root_to_group.into_values().collect()
    }
}