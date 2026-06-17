//! Embeddable DNA helix preview (character grid + RGB) for other TUIs.
//! The full-screen CLI remains in `main.rs`.

use std::f32::consts::TAU;
use std::time::Duration;

/// Target frame interval for helix mode in `main.rs` (`RenderState::draw` → helix branch).
const HELIX_FRAME_MS: f32 = 33.0;

pub type HelixCell = Option<(char, u8, u8, u8)>;
pub type HelixGrid = Vec<Vec<HelixCell>>;

#[allow(dead_code)]
const DNA: &[u8; 4] = b"ATCG";
const DNA_CODONS: &[&[u8; 3]] = &[
    b"ATG", b"GCT", b"GCC", b"GCA", b"GCG", b"TGT", b"TGC", b"GAT", b"GAC", b"GAA", b"GAG", b"TTT",
    b"TTC", b"GGT", b"GGC", b"GGA", b"GGG", b"CAT", b"CAC", b"ATT", b"ATC", b"ATA", b"AAA", b"AAG",
    b"TTA", b"TTG", b"CTT", b"CTC", b"CTA", b"CTG", b"AAT", b"AAC", b"CCT", b"CCC", b"CCA", b"CCG",
    b"CAA", b"CAG", b"CGT", b"CGC", b"CGA", b"CGG", b"AGA", b"AGG", b"TCT", b"TCC", b"TCA", b"TCG",
    b"AGT", b"AGC", b"ACT", b"ACC", b"ACA", b"ACG", b"GTT", b"GTC", b"GTA", b"GTG", b"TGG", b"TAT",
    b"TAC",
];

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn range(&mut self, max: u32) -> u32 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((self.0 >> 32) as u32) % max.max(1)
    }
}

fn put_cell(
    row: &mut [Option<(char, u8, u8, u8)>],
    width: usize,
    gx: i32,
    ch: u8,
    rgb: (u8, u8, u8),
) {
    if gx >= 1 && gx <= width as i32 {
        let ux = (gx - 1) as usize;
        if ux < width {
            row[ux] = Some((ch as char, rgb.0, rgb.1, rgb.2));
        }
    }
}

fn complement_dna(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base,
    }
}

fn coding_sequence_dna(codon_count: usize, rng: &mut Rng) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(codon_count * 3);
    sequence.extend_from_slice(b"ATG");
    for _ in 1..codon_count {
        sequence.extend_from_slice(DNA_CODONS[rng.range(DNA_CODONS.len() as u32) as usize]);
    }
    sequence
}

/// Classic-palette style colors (approximate ANSI mapping).
fn palette_front() -> (u8, u8, u8) {
    (0, 255, 255)
}

fn palette_back() -> (u8, u8, u8) {
    (0, 100, 255)
}

/// Stateful mini helix for embedding in a fixed character rectangle.
pub struct HelixMini {
    theta: f32,
    tick: usize,
    sequence: Vec<u8>,
}

impl HelixMini {
    pub fn new(seed: u64) -> Self {
        let mut rng = Rng::new(seed ^ 0xC0D0_5EED_B10C_0D0A);
        Self {
            theta: 0.0,
            tick: 0,
            // Same codon count as `HelixState` in `main.rs` for matching scroll length.
            sequence: coding_sequence_dna(720, &mut rng),
        }
    }

    /// Advance by one nominal helix frame (~33 ms at `speed` 1.0), matching `main.rs` timing.
    pub fn tick(&mut self, speed: f32) {
        self.tick_elapsed(speed, Duration::from_secs_f32(HELIX_FRAME_MS / 1000.0));
    }

    /// Advance using wall-clock time between draws so motion matches the CLI when the host
    /// UI redraws slowly. Same model as `main.rs`: per frame `theta += 0.16 * speed`,
    /// `tick += speed.ceil()`, frame sleep `33ms / speed`.
    pub fn tick_elapsed(&mut self, speed: f32, elapsed: Duration) {
        if speed <= 0.0 {
            return;
        }
        let elapsed_sec = elapsed.as_secs_f32().min(0.25);
        let frame_sec = (HELIX_FRAME_MS / 1000.0) / speed.max(0.01);
        let num_frames = elapsed_sec / frame_sec;
        self.theta += 0.16 * speed * num_frames;
        let tick_step = (speed.ceil() * num_frames).floor() as usize;
        self.tick = self.tick.wrapping_add(tick_step);
    }

    /// `grid[y][x]` = optional (base char, rgb). Unoccupied cells are `None`.
    pub fn render_grid(&self, width: usize, height: usize, scale: f32) -> HelixGrid {
        let mut grid: HelixGrid = vec![vec![None; width]; height];
        if width == 0 || height == 0 {
            return grid;
        }

        let mid_x = (width / 2).max(1) as i32;
        let rows = height.max(1);
        // Match `HelixState::draw` in `main.rs` (amplitude in terminal cells).
        let amp = (((width as f32) * 0.22) * scale).clamp(4.0, 48.0);

        for (yi, y) in (1_i32..=rows as i32).enumerate() {
            if yi >= height {
                break;
            }
            let t = y as f32 * 0.38 + self.theta;
            let z = t.sin();
            let x1 = mid_x + (t.cos() * amp) as i32;
            let x2 = mid_x + ((t + TAU / 2.0).cos() * amp) as i32;
            let base = self.sequence[(self.tick / 4 + y as usize) % self.sequence.len()];
            let pair = complement_dna(base);
            let (front_x, back_x, front_base, back_base) = if z >= 0.0 {
                (x1, x2, base, pair)
            } else {
                (x2, x1, pair, base)
            };

            let row = &mut grid[yi];

            let lo = back_x.min(front_x) + 1;
            let hi = back_x.max(front_x);
            for x in lo..hi {
                if x >= 1 && x <= width as i32 {
                    let ux = (x - 1) as usize;
                    if ux < width {
                        row[ux] = Some(('-', 90, 90, 90));
                    }
                }
            }

            put_cell(row, width, back_x, back_base, palette_back());
            put_cell(row, width, front_x, front_base, palette_front());
        }

        grid
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helix_grid_has_reasonable_extent() {
        let h = HelixMini::new(42);
        let g = h.render_grid(40, 10, 1.0);
        assert_eq!(g.len(), 10);
        assert!(g[0].len() == 40);
        let filled = g.iter().flatten().filter(|c| c.is_some()).count();
        assert!(filled > 5, "expected some helix cells");
    }
}
