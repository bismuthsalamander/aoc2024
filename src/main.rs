
use std::{collections::{HashMap, HashSet}, fs::read_to_string};

fn load_file_flat(filename: &str) -> String {
    read_to_string(filename).unwrap()
}

fn load_file_lines(filename: &str) -> Vec::<String> {
    read_to_string(filename).unwrap().lines().map(|s| String::from(s)).collect()
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

struct CharGrid {
    v: Vec::<char>,
    rows: usize,
    cols: usize
}

impl CharGrid {
    fn valid(&self, x: i64, y: i64) -> bool {
        x < self.cols as i64 && y < self.rows as i64 && x >= 0 && y >= 0
    }
    fn find(&self, c: char) -> (usize, usize) {
        for x in 0..self.cols {
            for y in 0..self.rows {
                if self.get(x as i64, y as i64) == c {
                    return (x, y);
                }
            }
        }
        panic!("Couldn't find {}", c);
    }
    fn get(&self, c: i64, y: i64) -> char {
        self.v[(y as usize*self.cols) + c as usize]
    }
    fn set(&mut self, x: i64, y: i64, c: char) {
        self.v[(y as usize*self.cols) + x as usize] = c;
    }
    fn check_match(&self, probe: &String, x: i64, y: i64, dx: i64, dy: i64) -> bool {
        for (offset, ch) in probe.chars().enumerate() {
            let oi64 = offset as i64;
            let x = x + (dx * oi64);
            let y = y + (dy * oi64);
            if !self.valid(x, y) {
                return false;
            }
            if self.get(x, y) != ch {
                return false;
            }
        }
        true
    }
}

fn load_grid(filename: &str) -> CharGrid {
    let mut v: Vec::<char> = Vec::new();
    let mut rows = 0;
    let mut cols = 0;
    for line in read_to_string(filename).unwrap().lines() {
        for char in line.chars() {
            if rows == 0 {
                cols += 1;
            }
            v.push(char);
        }
        rows += 1;
    }
    CharGrid {
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
    day4();
    day5();
    // day6();
    // day7();
    // day8();
    day9();
    day10();
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

fn day4() {
    println!("========Day 4========");
    let g = load_grid("day4.txt");
    let probe = String::from("XMAS");
    let mut ct = 0;
    for y in 0..g.rows {
        for x in 0..g.cols {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    if g.check_match(&probe, x as i64, y as i64, dx, dy) {
                        ct += 1;
                    }
                }
            }
        }
    }
    println!("XMAS count: {}", ct);

    let mut ct = 0;
    let probe = String::from("MAS");
    for y in 1i64..(g.rows - 1) as i64 {
        for x in 1i64..(g.cols - 1) as i64 {
            let mut found_probes = 0;
            for (dx, dy) in [(1, 1i64), (-1i64, 1)] {
                if g.check_match(&probe, x-dx, y-dy, dx, dy) {
                    found_probes += 1;
                } else if g.check_match(&probe, x+dx, y+dy, -dx, -dy) {
                    found_probes += 1;
                }
            }
            if found_probes >= 2 {
                ct += 1;
            }
        }
    }
    println!("X-MAS count: {}", ct);
}

fn day5() {
    println!("========Day 5========");
    let mut rules = HashSet::<(usize, usize)>::new();
    let lines = load_file_lines("day5.txt");
    let mut sum = 0;
    let mut sum_wrong = 0;
    let mut rules_done = false;
    for line in lines.iter() {
        if !rules_done {
            if line.len() == 0 {
                rules_done = true;
                continue;
            }
            let idx = line.find('|').unwrap();
            let first: usize = line[0..idx].parse().unwrap();
            let second: usize = line[idx+1..line.len()].parse().unwrap();
            rules.insert((first, second));
        } else {
            let mut pages: Vec<usize> = line.split(',').map(|s| s.parse::<usize>().unwrap()).collect();
            let mut ok = true;
            'outer: for second in 1..pages.len() {
                for first in 0..second {
                    if rules.contains(&(pages[second], pages[first])) {
                        ok = false;
                        break 'outer;
                    }
                }
            }
            if ok {
                sum += pages[(pages.len() - 1) / 2];
            } else {
                // Bubble sort - I'm not sure I could prove why we can get away with only examining
                // immediately adjacent elements, but it seems like a consequence of the ordering
                // being unique
                let mut changed = true;
                while changed {
                    changed = false;
                    for idx in 1..pages.len() {
                        if rules.contains(&(pages[idx], pages[idx-1])) {
                            pages.swap(idx, idx-1);
                            changed = true;
                        }
                    }
                }
                sum_wrong += pages[(pages.len() - 1) / 2];
            }
        }
    }
    println!("Sum of correctly ordered pages: {}", sum);
    println!("Sum of incorrectly ordered pages: {}", sum_wrong);
}

fn day6_run(g: &CharGrid) -> (usize, bool) {
    let mut visited = HashSet::<(i64, i64)>::new();
    let mut prev_states = HashSet::<(i64, i64, i64, i64)>::new();
    let (x , y) = g.find('^');
    let mut x = x as i64;
    let mut y = y as i64;
    let (mut dx, mut dy) = (0i64, -1i64);
    visited.insert((x, y));
    loop {
        if !g.valid(x + dx, y + dy) {
            return (visited.len(), false);
        }
        while g.get(x + dx, y + dy) == '#' {
            if prev_states.contains(&(x, y, dx, dy)) {
                return (visited.len(), true);
            }
            prev_states.insert((x, y, dx, dy));
            (dx, dy) = match (dx, dy) {
                (0, 1) => (-1, 0),
                (-1, 0) => (0, -1),
                (0, -1) => (1, 0),
                (1, 0) => (0, 1),
                _ => panic!("unknown direction")
            };
        }
        x += dx;
        y += dy;
        if g.get(x, y) == '#' {
            println!("x: {}, y: {}", x, y);
            panic!("whoops - stepped onto an obstacle");
        }
        visited.insert((x, y));
    }
}

fn day6() {
    println!("========Day 6========");
    let mut g = load_grid("day6.txt");
    let (visited, _luup) = day6_run(&g);
    println!("Locations: {}", visited);
    let mut luups = 0;
    for x in 0..g.cols as i64 {
        println!("col {} / {}", x, g.cols);
        for y in 0..g.rows as i64 {
            if g.get(x, y) != '.' {
                continue;
            }
            g.set(x, y, '#');
            let (_, luup) = day6_run(&g);
            if luup {
                luups += 1;
            }
            g.set(x, y, '.');
        }
    }
    println!("Luup locations: {}", luups);
}

fn day7_combo_update(v: &mut Vec<usize>, max: usize) -> bool {
    v[0] += 1;
    let mut idx = 0;
    while v[idx] >= max {
        v[idx] = 0;
        idx += 1;
        if idx == v.len() {
            return false;
        }
        v[idx] += 1;
    }
    true
}

fn day7() {
    println!("========Day 7========");
    let lines = load_file_lines("day7.txt");
    let mut problems = Vec::<(usize, Vec::<usize>)>::new();
    for line in lines.iter() {
        let colon = line.find(':').unwrap();
        let first: usize = line[0..colon].parse().unwrap();
        let p = (first, line[colon+2..line.len()].split(' ').map(|s| s.parse::<usize>().unwrap()).collect());
        problems.push(p);
    }
    let mut sum = 0;
    for p in problems.iter() {
        let mut ok = false;
        let max_mask = 2usize.pow(p.1.len() as u32);
        for mask in 0..max_mask {
            let mut val = p.1[0];
            for i in 0..(p.1.len() - 1) {
                if (mask & (1 << i)) > 0 {
                    val *= p.1[i+1];
                } else {
                    val += p.1[i+1];
                }
            }
            if val == p.0 {
                ok = true;
                break;
            }
        }
        if ok {
            sum += p.0;
        }
    }
    println!("Sum: {}", sum);
    let mut sum = 0;
    let mut z = 0;
    for p in problems.iter() {
        println!("{} {}", z, p.1.len());
        z += 1;
        let mut opers: Vec<usize> = (0..p.1.len()).map(|_idx| 0).collect();
        let mut ok = false;
        loop {
            let mut val = p.1[0];
            for i in 0..(p.1.len() - 1) {
                if opers[i] == 0 {
                    val += p.1[i+1];
                } else if opers[i] == 1 {
                    val *= p.1[i+1];
                } else if opers[i] == 2 {
                    let mut tmp = p.1[i+1];
                    while tmp > 0 {
                        val *= 10;
                        tmp /= 10;
                    }
                    val += p.1[i+1];
                } else {
                    panic!("whoops - val is {}", opers[i]);
                }
                if val > p.0 {
                    break;
                }
            }
            if val == p.0 {
                ok = true;
                break;
            }
            if !day7_combo_update(&mut opers, 3) {
                break;
            }
        }
        if ok {
            sum += p.0;
        }
    }
    println!("Sum: {}", sum);
}

fn day8() {
    println!("========Day 8========");
    let g = load_grid("day8.txt");   
    let mut antennas = HashMap::<char, Vec<(i64, i64)>>::new();
    for x in 0..g.cols as i64 {
        for y in 0..g.rows as i64 {
            let ch = g.get(x, y);
            if ch == '.' {
                continue;
            }
            antennas.entry(ch).or_insert(Vec::new()).push((x, y));
        }
    }

    let mut locations = HashSet::<(i64, i64)>::new();
    for (_k, coords) in antennas.iter() {
        for i in 0..coords.len() {
            for j in (i+1)..coords.len() {
                let (xi, yi) = coords[i];
                let (xj, yj) = coords[j];
                let dx = xj - xi;
                let dy = yj - yi;
                let x0 = xi - dx;
                let y0 = yi - dy;
                let x1 = xj + dx;
                let y1 = yj + dy;
                if g.valid(x0, y0) {
                    locations.insert((x0, y0));
                }
                if g.valid(x1, y1) {
                    locations.insert((x1, y1));
                }
            }
        }
    }
    println!("Locations: {}", locations.len());

    locations.clear();
    for (_k, coords) in antennas.iter() {
        for i in 0..coords.len() {
            for j in (i+1)..coords.len() {
                let (xi, yi) = coords[i];
                let (xj, yj) = coords[j];
                let dx = xj - xi;
                let dy = yj - yi;
                // we would need to divide dx and dy by their GCD, but each of the pairs in the
                // input data have dx and dy that are relatively prime.
                let mut x0 = xi;
                let mut y0 = yi;
                while g.valid(x0, y0) {
                    locations.insert((x0, y0));
                    x0 -= dx;
                    y0 -= dy;
                }
                x0 = xi;
                y0 = yi;
                while g.valid(x0, y0) {
                    locations.insert((x0, y0));
                    x0 += dx;
                    y0 += dy;
                }
            }
        }
    }
    // for y in 0i64..g.rows as i64 {
    //     for x in 0i64..g.cols as i64 {
    //         print!("{}", if locations.contains(&(x, y)) { '#' } else { g.get(x, y) });
    //     }
    //     println!("");
    // }
    println!("Locations: {}", locations.len());
}

fn day9() {
    println!("========Day 9========");
    let g = load_grid("day9.txt");
    let nums: Vec::<usize> = g.v.iter().map(|s| String::from(*s).parse().unwrap()).collect();
    let size: usize = nums.iter().sum();
    // None will represent a free block, and Some(n) represents a block for file n.
    let mut fs = Vec::<Option<usize>>::with_capacity(size);
    let mut is_file = true;
    let mut file_id = 0;
    for num in nums.iter() {
        for _i in 0..*num {
            fs.push(if is_file { Some(file_id) } else { None });
        }
        if is_file {
            file_id += 1;
        }
        is_file = !is_file;
    }
    let mut free_ptr = fs.iter().position(|val| val.is_none()).unwrap();
    let mut file_ptr = fs.iter().rposition(|val| val.is_some()).unwrap();
    while file_ptr > free_ptr {
        fs[free_ptr] = fs[file_ptr];
        fs[file_ptr] = None;
        while fs[file_ptr].is_none() {
            file_ptr -= 1;
        }
        while fs[free_ptr].is_some() {
            free_ptr += 1;
        }
    }
    // for num in fs.iter() {
    //     if num.is_none() {
    //         print!(".");
    //     } else {
    //         print!("{}", num.unwrap());
    //     }
    // }
    print!("\n");
    println!("Checksum: {}", day9_checksum(&fs));

    // Part 2 ********************************************
    // Okay, so I implemented it differently for part 2. I was annoyed that each search for a free
    // spot for a *single* file would be O(n^2) - I think, at least - in the original
    // implementation, so I did it this way. Fortunately, we don't have to worry about merging
    // free chunks since the moves all happen right-to-left.
    struct FSNode {
        size: usize,
        file_id: Option<usize>
    }
    impl FSNode {
        fn is_file(&self) -> bool {
            self.file_id.is_some()
        }
        fn is_free(&self) -> bool {
            !self.is_file()
        }
    }
    let mut fs = Vec::<FSNode>::with_capacity(nums.len());
    let mut is_file = true;
    let mut file_id = 0;
    for num in nums.iter() {
        fs.push(FSNode {
            file_id: if is_file { Some(file_id) } else { None },
            size: *num
        });
        if is_file {
            file_id += 1;
        }
        is_file = !is_file;
    }
    
    let mut first_free_ptr = fs.iter().position(|val| val.is_free()).unwrap();
    let mut file_ptr = fs.iter().rposition(|val| val.is_file()).unwrap();
    // We carry over the file_id variable to track the current search target, which prevents us
    // from mistakenly moving a file a second time. That wasn't a problem in part 1 because the
    // occupied space is completely contiguous, and we can abort the defrag as soon as the free
    // and file pointers cross over each other. Here, we may well encounter the same file twice
    // long before first_free_ptr and file_ptr meet.
    while first_free_ptr < file_ptr {
        let mut free_ptr = first_free_ptr;
        let end = fs[file_ptr].file_id == Some(0);
        while free_ptr < file_ptr && !(fs[free_ptr].is_free() && fs[free_ptr].size >= fs[file_ptr].size) {
            free_ptr += 1;
        }
        if fs[free_ptr].is_free() && fs[free_ptr].size >= fs[file_ptr].size {
            let leftover = fs[free_ptr].size - fs[file_ptr].size;
            fs[free_ptr].size = fs[file_ptr].size;
            fs[free_ptr].file_id = fs[file_ptr].file_id;
            fs[file_ptr].file_id = None;
            if leftover > 0 {
                fs.insert(free_ptr + 1, FSNode {
                    size: leftover,
                    file_id: None
                });
                file_ptr += 1;
            }
        }
        if end {
            break;
        }
        file_id -= 1;
        // find next file
        while file_ptr > 0 && fs[file_ptr].is_free() || fs[file_ptr].file_id.unwrap() >= file_id {
            file_ptr -= 1;
        }
        // find next free
        while fs[first_free_ptr].is_file() {
            first_free_ptr += 1;
        }
    }
    // for node in fs.iter() {
    //     for _ in 0..node.size {
    //         if node.is_file() {
    //             print!("{}", node.file_id.unwrap());
    //         } else {
    //             print!(".");
    //         }
    //     }
    // }
    let mut sum = 0;
    let mut block_id = 0;
    for node in fs.iter() {
        if node.is_free() {
            block_id += node.size;
        } else {
            for _ in 0..node.size {
                sum += block_id * node.file_id.unwrap();
                block_id += 1;
            }
        }
    }
    println!("Checksum: {}", sum);
}

fn day9_checksum(v: &Vec<Option<usize>>) -> usize {
    v.iter().enumerate().map(|(idx, x)| x.unwrap_or(0) * idx).sum()
}

fn day10_eval(g: &CharGrid, tx: i64, ty: i64) -> (usize, usize) {
    let mut dests = HashSet::<(i64, i64)>::new();
    let mut rating = 0;
    let mut frontier = Vec::<(i64, i64)>::new();
    frontier.push((tx, ty));
    while frontier.len() > 0 {
        let (x, y) = frontier.remove(0);
        let ch = g.get(x, y);
        if ch == '9' {
            dests.insert((x, y));
            // Luckily, we know that each visit to a given cell represents a unique path, so we
            // don't need to deduplicate anything for trailhead ratings.
            rating += 1;
        } else {
            for (dx, dy) in [(0, 1), (0, -1), (1, 0), (-1, 0)] {
                if g.valid(x+dx, y+dy) && g.get(x+dx, y+dy) == ((ch as u8) + 1) as char {
                    frontier.push((x+dx, y+dy));
                }
            }
        }
    }
    (dests.len(), rating)
}

fn day10() {
    println!("========Day 10========");
    let g = load_grid("day10.txt");
    let mut score = 0;
    let mut rating = 0;
    for y in 0..g.rows as i64 {
        for x in 0..g.cols as i64 {
            if g.get(x, y) != '0' {
                continue;
            }
            let (s, r) = day10_eval(&g, x, y);
            println!("{},{}", s, r);
            score += s;
            rating += r;
        }
    }
    println!("Score: {}\nRating: {}", score, rating);
}