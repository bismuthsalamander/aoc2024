
use std::{collections::HashMap, fs::read_to_string};

fn load_file_flat(filename: &str) -> String {
    read_to_string(filename).unwrap()
}

fn load_file(filename: &str) -> [Vec<i64>; 2] {
    let mut out = [Vec::new(), Vec::new()];
    for line in read_to_string(filename).unwrap().lines() {
        let mut idx = 0usize;
        for num in line.split(" ") {
            if num.len() == 0 {
                continue;
            }
            out[idx].push(String::from(num).parse().unwrap());
            idx += 1;
        }
    }
    out[0].sort();
    out[1].sort();
    out
}

pub fn day1() {
    println!("========Day 1========");
    let locations = load_file("day1.txt");
    let mut dist = 0;
    for idx in 0..locations[0].len() {
        dist += (locations[0][idx] - locations[1][idx]).abs();
    }
    println!("Dist: {}", dist);
    let mut right_counts = HashMap::new();
    for num in locations[1].iter() {
        *right_counts.entry(num).or_insert(0) += 1;
    }
    let mut similarity = 0;
    for num in locations[0].iter() {
        similarity += num * right_counts.get(&num).or(Some(&0)).unwrap();
    }
    println!("Similarity: {}", similarity);
}

struct Grid {
    v: Vec::<i64>,
    rows: usize,
    cols: usize
}

impl Grid {
    fn get(&self, y: usize, c: usize) -> i64 {
        self.v[y*self.cols + c]
    }
}

fn load_grid(filename: &str) -> Grid {
    let mut v: Vec::<i64> = Vec::new();
    let mut rows = 0;
    let mut cols = 0;
    for line in read_to_string(filename).unwrap().lines() {
        for num in line.split(" ") {
            if num.len() == 0 {
                continue;
            }
            if rows == 0 {
                cols += 1;
            }
            v.push(String::from(num).parse().unwrap());
        }
        rows += 1;
    }
    Grid {
        v,
        rows,
        cols
    }
}

fn load_vec(filename: &str) -> Vec<Vec<i64>> {
    let mut v: Vec::<Vec::<i64>> = Vec::new();
    for line in read_to_string(filename).unwrap().lines() {
        let mut row = Vec::new();
        for num in line.split(" ") {
            if num.len() == 0 {
                continue;
            }
            row.push(String::from(num).parse().unwrap());
        }
        v.push(row);
    }
    v
}

fn day2() {
    println!("========Day 2========");
    let v = load_vec("day2.txt");
    day2_pt1(&v);
    day2_pt2(&v);
}

fn day2_pt1(v: &Vec<Vec<i64>>) {
    let mut safe_ct = 0;
    for row in v.iter() {
        let (ok, _) = day2_safe(&row);
        if ok {
            safe_ct += 1;
        }
    }
    println!("Safe: {}", safe_ct);
}

fn day2_safe(row: &Vec<i64>) -> (bool, usize) {
    if row.len() < 2 {
        return (true, 0);
    }    
    let inc = row[1] > row[0];
    for c in 1..row.len() {
        let diff = row[c] - row[c-1];
        if (diff > 0) != inc {
            return (false, c);
        }
        if diff.abs() > 3 || diff.abs() < 1 {
            return (false, c);
        }
    }
    return (true, 0);
}

fn day2_safe_skip(row: &Vec<i64>, skip: usize) -> (bool, usize) {
    if row.len() < 2 {
        return (true, 0);
    }
    let inc = row[if skip <= 1 { 2 } else { 1 } ] > row[if skip == 0 { 1 } else { 0 } ];
    
    let start = 1 + if skip < 2 { 1 } else { 0 };
    let mut offset = if start == skip + 1 { 1 } else { 0 };
    for c in start..row.len() {
        if c == skip {
            offset += 1;
            continue;
        }
        let diff = row[c] - row[c-(1 + offset)];
        if offset == 1 {
            offset -= 1;
        }
        if (diff > 0) != inc || diff.abs() > 3 || diff.abs() < 1 {
            return (false, c);
        }
        
    }
    (true, 0)
}

fn day2_pt2(v: &Vec<Vec<i64>>) {
    let mut safe_ct = 0;
    let mut _orig_safe_ct = 0;
    for (_ri, row )in v.iter().enumerate() {
        let (ok, _idx) = day2_safe(&row);
        if ok {
            safe_ct += 1;
            _orig_safe_ct += 1;
            continue;
        }
        let mut skip_res = false;
        let mut _skip_idx = 0;
        for skip in 0..row.len() {
            let mut new_row = row.clone();
            new_row.remove(skip);
            let (ok, _idx) = day2_safe(&new_row);
            if ok {
                skip_res = true;
                _skip_idx = skip;
                break;
            }
        }
        if skip_res {
            safe_ct += 1;
        }
    }
    println!("Safe with dampener: {}", safe_ct);
}

fn main() {
    day1();
    day2();
    day3();
}

fn day3_try_probe(s: String) -> Option<i32> {
    let mut state = 0;
    let mut offset = 0;
    let mut nums = [0, 0];
    for ch in s.chars() {
        if ch == ',' {
            if offset == 0 || offset > 3 {
                return None;
            }
            if state >= 1 {
                return None;
            }
            offset = 0;
            state = 1;
        } else if ch.is_numeric() {
            offset += 1;
            if offset > 3 {
                return None;
            }
            nums[state] *= 10;
            nums[state] += (ch as i32 - '0' as i32);
        } else if ch == ')' {
            return Some(nums[0] * nums[1]);
        } else {
            return None;
        }
    }
    return None;
}

fn day3() {
    println!("========Day 3========");
    let s = load_file_flat("day3.txt");
    let mut out = 0;
    for probe in s.split("mul(") {
        if let Some(product) = day3_try_probe(String::from(probe)) {
            out += product;
        }
    }
    println!("Sum: {}", out);

    let mut out = 0;
    for (idx, dont_pair) in s.split("don't()").enumerate() {
        let probes = if idx == 0 {
            dont_pair.split("mul(")
        } else {
            let Some((_before, after)) = dont_pair.split_once("do()") else {
                continue;
            };
            after.split("mul(")
        };
                
        for probe in probes {
            if let Some(product) = day3_try_probe(String::from(probe)) {
                out += product;
            }
        }
    }
    println!("Sum with do/don't: {}", out);
}