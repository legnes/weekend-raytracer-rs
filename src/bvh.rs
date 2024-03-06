use super::hit::HitRecord;
use super::hit::World;
use super::ray::Ray;
use super::vec::Point3;

fn ray_intersects_aabb(
    aabb_min: Point3,
    aabb_max: Point3,
    ray: &Ray,
    t_min: f64,
    t_max: f64,
) -> bool {
    let t1 = (aabb_min - ray.origin()) / ray.direction();
    let t2 = (aabb_max - ray.origin()) / ray.direction();
    let aabb_t_min = (t1.x().min(t2.x()))
        .max(t1.y().min(t2.y()))
        .max(t1.z().min(t2.z()));
    let aabb_t_max = (t1.x().max(t2.x()))
        .min(t1.y().max(t2.y()))
        .min(t1.z().max(t2.z()));
    aabb_t_max >= aabb_t_min && aabb_t_min < t_max && aabb_t_max > t_min
}

pub struct Bvh {
    // SE TODO: Can we remove this?
    root_index: usize,
    nodes_used: usize,
    // NOTE: This is actually fixed size (size of primitives ** 2 - 1)?
    nodes: Vec<Node>,
    // NOTE: This is also fixed size (size of the primitives)
    primitive_map: Vec<usize>,
}

impl Bvh {
    pub fn new(primitives: &World) -> Self {
        let n = primitives.len();
        let max_nodes = n.pow(2) - 1;

        let mut bvh = Self {
            root_index: 0,
            nodes_used: 1,
            nodes: std::iter::repeat_with(|| Node::new())
                .take(max_nodes)
                .collect::<Vec<_>>(),
            primitive_map: (0..n).collect(),
        };

        let root = &mut bvh.nodes[bvh.root_index];
        root.set_primitive_range(0, n);
        root.update_bounds(&bvh.primitive_map, primitives);

        bvh.subdivide_node(bvh.root_index, primitives);

        eprintln!("built bvh with {} nodes", bvh.nodes_used);
        bvh
    }

    fn subdivide_node(&mut self, node_index: usize, primitives: &World) -> () {
        // Split along longest axis
        let node = &mut self.nodes[node_index];
        let extent = node.aabb_max - node.aabb_min;
        let mut axis: usize = 0;
        if extent[1] > extent[0] {
            axis = 1;
        }
        if extent[2] > extent[axis] {
            axis = 2;
        }
        let split_position = node.aabb_min[axis] + extent[axis] * 0.5;

        // Bucket the primitive map into left and right
        // so that each node can reference a continuous slice of primitive indices
        let mut left_cursor = node.primitive_start;
        let mut right_cursor = left_cursor + node.primitive_count - 1;
        while left_cursor < right_cursor {
            if primitives[self.primitive_map[left_cursor]].centroid()[axis] < split_position {
                left_cursor += 1;
            } else {
                // Swap
                let temp = self.primitive_map[left_cursor];
                self.primitive_map[left_cursor] = self.primitive_map[right_cursor];
                self.primitive_map[right_cursor] = temp;
                right_cursor -= 1;
            }
        }

        let left_primitive_count = left_cursor - node.primitive_start;
        // Base case: only one node (all left or all right)
        if left_primitive_count == 0 || left_primitive_count == node.primitive_count {
            return ();
        }

        // Create child nodes
        node.left_child_index = self.nodes_used;
        self.nodes_used += 2;
        let primitive_count = node.primitive_count;
        node.primitive_count = 0;

        // Redeclare variables so node can go out of scope
        let left_child_index = node.left_child_index();
        let right_child_index = node.right_child_index();
        let primitive_start = node.primitive_start;

        let left_child = &mut self.nodes[left_child_index];
        left_child.set_primitive_range(primitive_start, left_primitive_count);
        left_child.update_bounds(&self.primitive_map, primitives);
        self.subdivide_node(left_child_index, primitives);

        let right_child = &mut self.nodes[right_child_index];
        right_child.set_primitive_range(left_cursor, primitive_count - left_primitive_count);
        right_child.update_bounds(&self.primitive_map, primitives);
        self.subdivide_node(right_child_index, primitives);
    }

    pub fn hit(&self, ray: &Ray, t_min: f64, t_max: f64, primitives: &World) -> Option<HitRecord> {
        self.hit_node(ray, t_min, t_max, primitives, self.root_index)
    }

    fn hit_node(
        &self,
        ray: &Ray,
        t_min: f64,
        t_max: f64,
        primitives: &World,
        node_index: usize,
    ) -> Option<HitRecord> {
        let node = &self.nodes[node_index];
        let mut closest_t = t_max;
        let mut closest_record = None;

        if !ray_intersects_aabb(node.aabb_min, node.aabb_max, ray, t_min, t_max) {
            return None;
        }

        if node.is_leaf() {
            // Hit all primitives on the node
            for i in node.primitive_start..(node.primitive_start + node.primitive_count) {
                let primitive = &primitives[self.primitive_map[i]];
                if let Some(record) = primitive.hit(ray, t_min, closest_t) {
                    closest_t = record.t;
                    closest_record = Some(record);
                }
            }
        } else {
            // Hit left node
            if let Some(record) =
                self.hit_node(ray, t_min, closest_t, primitives, node.left_child_index())
            {
                closest_t = record.t;
                closest_record = Some(record);
            }
            // Hit right node
            if let Some(record) =
                self.hit_node(ray, t_min, closest_t, primitives, node.right_child_index())
            {
                closest_record = Some(record);
            }
        }

        closest_record
    }
}

#[derive(Debug)]
pub struct Node {
    aabb_min: Point3,
    aabb_max: Point3,
    left_child_index: usize,
    primitive_start: usize,
    primitive_count: usize,
}

impl Node {
    pub fn new() -> Self {
        Self {
            aabb_min: Point3::zero(),
            aabb_max: Point3::zero(),
            left_child_index: 0,
            primitive_start: 0,
            primitive_count: 0,
        }
    }

    pub fn left_child_index(&self) -> usize {
        self.left_child_index
    }

    pub fn right_child_index(&self) -> usize {
        self.left_child_index + 1
    }

    pub fn is_leaf(&self) -> bool {
        self.primitive_count > 0
    }

    pub fn set_primitive_range(&mut self, start: usize, count: usize) -> () {
        self.primitive_start = start;
        self.primitive_count = count;
    }

    pub fn update_bounds(&mut self, primitive_map: &Vec<usize>, primitives: &World) -> () {
        self.aabb_min = Point3::fill(f64::INFINITY);
        self.aabb_max = Point3::fill(f64::NEG_INFINITY);
        for index_offset in 0..self.primitive_count {
            let primitive_index = primitive_map[self.primitive_start + index_offset];
            let primitive = &primitives[primitive_index];
            self.aabb_min = self.aabb_min.min(primitive.aabb_min());
            self.aabb_max = self.aabb_max.max(primitive.aabb_max());
        }
    }
}
