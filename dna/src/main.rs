use std::{
    env,
    f32::consts::TAU,
    io::{self, IsTerminal, Read, Write},
    process::Command,
    sync::mpsc::{self, Receiver},
    thread,
    time::{Duration, Instant},
};

const DNA: &[u8; 4] = b"ATCG";
const RNA: &[u8; 4] = b"AUCG";
const DNA_CODONS: &[&[u8; 3]] = &[
    b"ATG", b"GCT", b"GCC", b"GCA", b"GCG", b"TGT", b"TGC", b"GAT", b"GAC", b"GAA", b"GAG", b"TTT",
    b"TTC", b"GGT", b"GGC", b"GGA", b"GGG", b"CAT", b"CAC", b"ATT", b"ATC", b"ATA", b"AAA", b"AAG",
    b"TTA", b"TTG", b"CTT", b"CTC", b"CTA", b"CTG", b"AAT", b"AAC", b"CCT", b"CCC", b"CCA", b"CCG",
    b"CAA", b"CAG", b"CGT", b"CGC", b"CGA", b"CGG", b"AGA", b"AGG", b"TCT", b"TCC", b"TCA", b"TCG",
    b"AGT", b"AGC", b"ACT", b"ACC", b"ACA", b"ACG", b"GTT", b"GTC", b"GTA", b"GTG", b"TGG", b"TAT",
    b"TAC",
];
const RNA_CODONS: &[&[u8; 3]] = &[
    b"AUG", b"GCU", b"GCC", b"GCA", b"GCG", b"UGU", b"UGC", b"GAU", b"GAC", b"GAA", b"GAG", b"UUU",
    b"UUC", b"GGU", b"GGC", b"GGA", b"GGG", b"CAU", b"CAC", b"AUU", b"AUC", b"AUA", b"AAA", b"AAG",
    b"UUA", b"UUG", b"CUU", b"CUC", b"CUA", b"CUG", b"AAU", b"AAC", b"CCU", b"CCC", b"CCA", b"CCG",
    b"CAA", b"CAG", b"CGU", b"CGC", b"CGA", b"CGG", b"AGA", b"AGG", b"UCU", b"UCC", b"UCA", b"UCG",
    b"AGU", b"AGC", b"ACU", b"ACC", b"ACA", b"ACG", b"GUU", b"GUC", b"GUA", b"GUG", b"UGG", b"UAU",
    b"UAC",
];

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    let cli = match Cli::parse(&args) {
        Ok(cli) => cli,
        Err(message) => {
            eprintln!("{message}\n");
            print_help();
            return Ok(());
        }
    };
    if cli.help {
        print_help();
        return Ok(());
    }

    match cli.mode {
        #[cfg(feature = "codon-wheel")]
        RequestedMode::Codon => run_interactive(AppMode::Codon),
        #[cfg(not(feature = "codon-wheel"))]
        RequestedMode::Codon => {
            eprintln!("codon wheel is disabled in the default build");
            eprintln!("run with: cargo run --bin DNA --features codon-wheel -- codon");
            Ok(())
        }
        RequestedMode::Matrix => run_interactive(if cli.alphabet == Alphabet::Rna {
            AppMode::MatrixRna
        } else {
            AppMode::MatrixDna
        }),
        RequestedMode::Helix => run_interactive(if cli.alphabet == Alphabet::Rna {
            AppMode::HelixRna
        } else {
            AppMode::HelixDna
        }),
    }
}

fn print_help() {
    eprintln!(
        "DNA\n\nUSAGE:\n  DNA [MODE] [ALPHABET]\n\nMODES:\n  h, helix       Rotating double helix. This is the default.\n  m, matrix      Nucleotide rain, cmatrix-style.\n  codon          Experimental codon wheel. Requires --features codon-wheel.\n\nALPHABET:\n  d, dna, -d, --dna      Use ATCG. This is the default.\n  r, rna, -r, --rna      Use AUCG.\n\nEXAMPLES:\n  DNA              # DNA helix\n  DNA h r          # RNA helix\n  DNA m d          # DNA matrix\n  DNA matrix -r    # RNA matrix\n  cargo run --bin DNA -- h --rna\n  cargo run --bin DNA --features codon-wheel -- codon\n\nCONTROLS:\n  Left/Right    Cycle mode while running\n  Up/Down       Change speed; Down reaches 0.00x freeze\n  +/-           Change visual scale\n  c             Cycle color palettes\n  f             Toggle focus mode and hide footer\n  q or Ctrl-C   Quit cleanly"
    );
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Alphabet {
    Dna,
    Rna,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RequestedMode {
    Helix,
    Matrix,
    Codon,
}

struct Cli {
    mode: RequestedMode,
    alphabet: Alphabet,
    help: bool,
}

impl Cli {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut mode = RequestedMode::Helix;
        let mut alphabet = Alphabet::Dna;
        let mut mode_set = false;

        for arg in args {
            match arg.as_str() {
                "-h" | "--help" | "help" => {
                    return Ok(Self {
                        mode,
                        alphabet,
                        help: true,
                    });
                }
                "h" | "helix" => {
                    if mode_set {
                        return Err(format!("multiple modes supplied; unexpected `{arg}`"));
                    }
                    mode = RequestedMode::Helix;
                    mode_set = true;
                }
                "m" | "matrix" | "rain" => {
                    if mode_set {
                        return Err(format!("multiple modes supplied; unexpected `{arg}`"));
                    }
                    mode = RequestedMode::Matrix;
                    mode_set = true;
                }
                "codon" | "codons" | "table" => {
                    if mode_set {
                        return Err(format!("multiple modes supplied; unexpected `{arg}`"));
                    }
                    mode = RequestedMode::Codon;
                    mode_set = true;
                }
                "d" | "dna" | "-d" | "--dna" => alphabet = Alphabet::Dna,
                "r" | "rna" | "-r" | "--rna" => alphabet = Alphabet::Rna,
                _ => return Err(format!("unknown argument `{arg}`")),
            }
        }

        Ok(Self {
            mode,
            alphabet,
            help: false,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum AppMode {
    MatrixDna,
    MatrixRna,
    HelixDna,
    HelixRna,
    #[cfg(feature = "codon-wheel")]
    Codon,
}

impl AppMode {
    fn next(self) -> Self {
        match self {
            Self::MatrixDna => Self::MatrixRna,
            Self::MatrixRna => Self::HelixDna,
            Self::HelixDna => Self::HelixRna,
            #[cfg(feature = "codon-wheel")]
            Self::HelixRna => Self::Codon,
            #[cfg(not(feature = "codon-wheel"))]
            Self::HelixRna => Self::MatrixDna,
            #[cfg(feature = "codon-wheel")]
            Self::Codon => Self::MatrixDna,
        }
    }

    fn prev(self) -> Self {
        match self {
            #[cfg(feature = "codon-wheel")]
            Self::MatrixDna => Self::Codon,
            #[cfg(not(feature = "codon-wheel"))]
            Self::MatrixDna => Self::HelixRna,
            Self::MatrixRna => Self::MatrixDna,
            Self::HelixDna => Self::MatrixRna,
            Self::HelixRna => Self::HelixDna,
            #[cfg(feature = "codon-wheel")]
            Self::Codon => Self::HelixRna,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::MatrixDna => "DNA matrix",
            Self::MatrixRna => "RNA matrix",
            Self::HelixDna => "DNA helix",
            Self::HelixRna => "RNA helix",
            #[cfg(feature = "codon-wheel")]
            Self::Codon => "RNA codon wheel",
        }
    }

    fn alphabet(self) -> &'static [u8; 4] {
        match self {
            #[cfg(feature = "codon-wheel")]
            Self::MatrixRna | Self::HelixRna | Self::Codon => RNA,
            #[cfg(not(feature = "codon-wheel"))]
            Self::MatrixRna | Self::HelixRna => RNA,
            Self::MatrixDna | Self::HelixDna => DNA,
        }
    }
}

enum Key {
    Up,
    Down,
    Left,
    Right,
    Bigger,
    Smaller,
    Focus,
    Color,
    Quit,
}

#[derive(Clone, Copy)]
enum Palette {
    Classic,
    Bases,
    Ice,
    Fire,
    Mono,
}

impl Palette {
    fn next(self) -> Self {
        match self {
            Self::Classic => Self::Bases,
            Self::Bases => Self::Ice,
            Self::Ice => Self::Fire,
            Self::Fire => Self::Mono,
            Self::Mono => Self::Classic,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Classic => "classic",
            Self::Bases => "base colors",
            Self::Ice => "ice",
            Self::Fire => "fire",
            Self::Mono => "mono",
        }
    }

    fn matrix(self, shade: usize) -> &'static str {
        match (self, shade) {
            (_, 0) => "\x1b[97m",
            (Self::Classic, 1) => "\x1b[92m",
            (Self::Classic, _) => "\x1b[32m",
            (Self::Bases, 1) => "\x1b[97m",
            (Self::Bases, _) => "\x1b[37m",
            (Self::Ice, 1) => "\x1b[96m",
            (Self::Ice, _) => "\x1b[36m",
            (Self::Fire, 1) => "\x1b[93m",
            (Self::Fire, _) => "\x1b[31m",
            (Self::Mono, 1) => "\x1b[37m",
            (Self::Mono, _) => "\x1b[90m",
        }
    }

    fn back(self) -> &'static str {
        match self {
            Self::Classic => "\x1b[34m",
            Self::Bases => "\x1b[90m",
            Self::Ice => "\x1b[36m",
            Self::Fire => "\x1b[31m",
            Self::Mono => "\x1b[90m",
        }
    }

    fn front(self) -> &'static str {
        match self {
            Self::Classic => "\x1b[96m",
            Self::Bases => "\x1b[97m",
            Self::Ice => "\x1b[97m",
            Self::Fire => "\x1b[93m",
            Self::Mono => "\x1b[37m",
        }
    }

    #[cfg(feature = "codon-wheel")]
    fn codon(self, ring: usize, active: bool) -> &'static str {
        if active {
            return match (self, ring) {
                (Self::Classic, 0) => "\x1b[30;102m",
                (Self::Classic, 1) => "\x1b[30;106m",
                (Self::Classic, _) => "\x1b[30;103m",
                (Self::Bases, 0) => "\x1b[30;101m",
                (Self::Bases, 1) => "\x1b[30;104m",
                (Self::Bases, _) => "\x1b[30;103m",
                (Self::Ice, 0) => "\x1b[30;106m",
                (Self::Ice, 1) => "\x1b[30;104m",
                (Self::Ice, _) => "\x1b[30;107m",
                (Self::Fire, 0) => "\x1b[30;103m",
                (Self::Fire, 1) => "\x1b[30;101m",
                (Self::Fire, _) => "\x1b[30;105m",
                (Self::Mono, _) => "\x1b[30;107m",
            };
        }
        match (self, ring) {
            (Self::Classic, 0) => "\x1b[92m",
            (Self::Classic, 1) => "\x1b[96m",
            (Self::Classic, _) => "\x1b[93m",
            (Self::Bases, 0) => "\x1b[91m",
            (Self::Bases, 1) => "\x1b[94m",
            (Self::Bases, _) => "\x1b[93m",
            (Self::Ice, 0) => "\x1b[36m",
            (Self::Ice, 1) => "\x1b[94m",
            (Self::Ice, _) => "\x1b[97m",
            (Self::Fire, 0) => "\x1b[93m",
            (Self::Fire, 1) => "\x1b[91m",
            (Self::Fire, _) => "\x1b[95m",
            (Self::Mono, 0) => "\x1b[37m",
            (Self::Mono, 1) => "\x1b[90m",
            (Self::Mono, _) => "\x1b[97m",
        }
    }

    #[cfg(feature = "codon-wheel")]
    fn accent(self) -> &'static str {
        match self {
            Self::Classic => "\x1b[91m",
            Self::Bases => "\x1b[97m",
            Self::Ice => "\x1b[96m",
            Self::Fire => "\x1b[93m",
            Self::Mono => "\x1b[97m",
        }
    }

    fn base(self, base: u8) -> &'static str {
        if !matches!(self, Self::Bases) {
            return self.front();
        }
        match base {
            b'A' => "\x1b[91m",
            b'C' => "\x1b[94m",
            b'G' => "\x1b[93m",
            b'T' | b'U' => "\x1b[92m",
            _ => "\x1b[97m",
        }
    }

    #[cfg(feature = "codon-wheel")]
    fn active_base(self, base: u8) -> &'static str {
        if !matches!(self, Self::Bases) {
            return self.accent();
        }
        match base {
            b'A' => "\x1b[30;101m",
            b'C' => "\x1b[30;104m",
            b'G' => "\x1b[30;103m",
            b'T' | b'U' => "\x1b[30;102m",
            _ => "\x1b[30;107m",
        }
    }
}

struct TerminalMode {
    raw: bool,
}

impl TerminalMode {
    fn enter() -> Self {
        let raw = io::stdin().is_terminal()
            && Command::new("stty")
                .args(["raw", "-echo"])
                .status()
                .is_ok_and(|status| status.success());
        Self { raw }
    }
}

impl Drop for TerminalMode {
    fn drop(&mut self) {
        if self.raw {
            let _ = Command::new("stty").args(["sane"]).status();
        }
        let _ = write!(io::stdout(), "\x1b[0m\x1b[?25h\x1b[2J\x1b[H");
    }
}

fn spawn_input_thread() -> Receiver<Key> {
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let mut stdin = io::stdin().lock();
        let mut byte = [0_u8; 1];
        loop {
            match stdin.read(&mut byte) {
                Ok(0) => thread::sleep(Duration::from_millis(20)),
                Ok(_) => match byte[0] {
                    b'q' | 3 => {
                        let _ = tx.send(Key::Quit);
                        break;
                    }
                    b'+' | b'=' => {
                        let _ = tx.send(Key::Bigger);
                    }
                    b'-' | b'_' => {
                        let _ = tx.send(Key::Smaller);
                    }
                    b'f' | b'F' => {
                        let _ = tx.send(Key::Focus);
                    }
                    b'c' | b'C' => {
                        let _ = tx.send(Key::Color);
                    }
                    27 => {
                        let mut seq = [0_u8; 2];
                        if stdin.read_exact(&mut seq).is_ok() && seq[0] == b'[' {
                            let key = match seq[1] {
                                b'A' => Some(Key::Up),
                                b'B' => Some(Key::Down),
                                b'C' => Some(Key::Right),
                                b'D' => Some(Key::Left),
                                _ => None,
                            };
                            if let Some(key) = key {
                                let _ = tx.send(key);
                            }
                        }
                    }
                    _ => {}
                },
                Err(_) => break,
            }
        }
    });
    rx
}

fn run_interactive(initial_mode: AppMode) -> io::Result<()> {
    let _terminal = TerminalMode::enter();
    let input = spawn_input_thread();
    let mut out = io::stdout().lock();
    let mut mode = initial_mode;
    let seed = 1_u64;
    let mut speed = 1.0_f32;
    let mut scale = 1.0_f32;
    let mut focus = false;
    let mut palette = Palette::Classic;
    let mut state = RenderState::new(mode, seed);

    write!(out, "\x1b[?25l")?;
    clear(&mut out)?;
    loop {
        let start = Instant::now();
        while let Ok(key) = input.try_recv() {
            match key {
                Key::Quit => return Ok(()),
                Key::Left => {
                    mode = mode.prev();
                    state = RenderState::new(mode, seed);
                    clear(&mut out)?;
                }
                Key::Right => {
                    mode = mode.next();
                    state = RenderState::new(mode, seed);
                    clear(&mut out)?;
                }
                Key::Up => {
                    speed = if speed == 0.0 {
                        0.15
                    } else {
                        (speed * 1.25).min(4.0)
                    };
                }
                Key::Down => {
                    speed = if speed <= 0.15 { 0.0 } else { speed / 1.25 };
                }
                Key::Bigger => {
                    scale = (scale * 1.15).min(2.5);
                }
                Key::Smaller => {
                    scale = (scale / 1.15).max(0.5);
                }
                Key::Focus => {
                    focus = !focus;
                    clear(&mut out)?;
                }
                Key::Color => {
                    palette = palette.next();
                    clear(&mut out)?;
                }
            }
        }

        let frame = state.draw(&mut out, mode, seed, speed, scale, focus, palette)?;
        out.flush()?;
        sleep_frame(start, scaled_frame(frame, speed));
    }
}

enum RenderState {
    Matrix(MatrixState),
    Helix(HelixState),
    #[cfg(feature = "codon-wheel")]
    Codon(CodonState),
}

impl RenderState {
    fn new(mode: AppMode, seed: u64) -> Self {
        match mode {
            AppMode::MatrixDna | AppMode::MatrixRna => {
                Self::Matrix(MatrixState::new(seed, mode.alphabet()))
            }
            AppMode::HelixDna | AppMode::HelixRna => {
                Self::Helix(HelixState::new(seed, mode.alphabet()))
            }
            #[cfg(feature = "codon-wheel")]
            AppMode::Codon => Self::Codon(CodonState::new(seed)),
        }
    }

    fn draw(
        &mut self,
        out: &mut impl Write,
        mode: AppMode,
        seed: u64,
        speed: f32,
        scale: f32,
        focus: bool,
        palette: Palette,
    ) -> io::Result<Duration> {
        match self {
            Self::Matrix(state) => {
                state.draw(out, mode, seed, speed, scale, focus, palette)?;
                Ok(Duration::from_millis(45))
            }
            Self::Helix(state) => {
                state.draw(out, mode, seed, speed, scale, focus, palette)?;
                Ok(Duration::from_millis(33))
            }
            #[cfg(feature = "codon-wheel")]
            Self::Codon(state) => {
                state.draw(out, seed, speed, scale, focus, palette)?;
                Ok(Duration::from_millis(320))
            }
        }
    }
}

struct MatrixState {
    width: u16,
    height: u16,
    rng: Rng,
    columns: Columns,
}

impl MatrixState {
    fn new(seed: u64, alphabet: &[u8; 4]) -> Self {
        let (width, height) = terminal_size();
        let mut rng = Rng::new(seed ^ alphabet_seed(alphabet) ^ 0xA7C9_51E5_DA7A_B105);
        let columns = Columns::new(width as usize, height as i32, &mut rng);
        Self {
            width,
            height,
            rng,
            columns,
        }
    }

    fn draw(
        &mut self,
        out: &mut impl Write,
        mode: AppMode,
        seed: u64,
        speed: f32,
        scale: f32,
        focus: bool,
        palette: Palette,
    ) -> io::Result<()> {
        let alphabet = mode.alphabet();
        let size = terminal_size();
        if size != (self.width, self.height) {
            (self.width, self.height) = size;
            self.columns = Columns::new(self.width as usize, self.height as i32, &mut self.rng);
        }

        let column_step = if scale >= 1.8 {
            3
        } else if scale >= 1.25 {
            2
        } else {
            1
        };

        if speed == 0.0 {
            return status_line(out, self.height, mode, seed, speed, scale, focus, palette);
        }

        write!(out, "\x1b[2J")?;
        for x in (0..self.width as usize).step_by(column_step) {
            let col = &mut self.columns.0[x];
            for i in 0..col.len {
                let y = col.head - i;
                if y < 1 || y > self.height as i32 {
                    continue;
                }
                let base = alphabet[self.rng.range(4) as usize];
                let shade = if matches!(palette, Palette::Bases) {
                    palette.base(base)
                } else if i == 0 {
                    palette.matrix(0)
                } else if i < 4 {
                    palette.matrix(1)
                } else {
                    palette.matrix(2)
                };
                write!(out, "\x1b[{y};{}H{shade}{}", x + 1, base as char)?;
            }
            col.head += col.speed;
            if col.head - col.len > self.height as i32 {
                col.reset(self.height as i32, &mut self.rng);
            }
        }
        status_line(out, self.height, mode, seed, speed, scale, focus, palette)
    }
}

struct HelixState {
    theta: f32,
    tick: usize,
    sequence: Vec<u8>,
}

impl HelixState {
    fn new(seed: u64, alphabet: &[u8; 4]) -> Self {
        let mut rng = Rng::new(seed ^ alphabet_seed(alphabet) ^ 0xC0D0_5EED_B10C_0D0A);
        Self {
            theta: 0.0,
            tick: 0,
            sequence: coding_sequence(alphabet, 720, &mut rng),
        }
    }

    fn draw(
        &mut self,
        out: &mut impl Write,
        mode: AppMode,
        seed: u64,
        speed: f32,
        scale: f32,
        focus: bool,
        palette: Palette,
    ) -> io::Result<()> {
        let alphabet = mode.alphabet();
        let (width, height) = terminal_size();
        let mid_x = (width / 2).max(1) as i32;
        let rows = height.saturating_sub(2).max(1) as i32;
        let amp = (((width as f32) * 0.22) * scale).clamp(4.0, 48.0);

        clear(out)?;
        for y in 1..=rows {
            let t = y as f32 * 0.38 + self.theta;
            let z = t.sin();
            let x1 = mid_x + (t.cos() * amp) as i32;
            let x2 = mid_x + ((t + TAU / 2.0).cos() * amp) as i32;
            let base = self.sequence[(self.tick / 4 + y as usize) % self.sequence.len()];
            let pair = complement(base, alphabet);
            let (front_x, back_x, front_base, back_base) = if z >= 0.0 {
                (x1, x2, base, pair)
            } else {
                (x2, x1, pair, base)
            };

            draw_bridge(out, y, back_x, front_x)?;
            draw_base(
                out,
                y,
                back_x,
                back_base as char,
                if matches!(palette, Palette::Bases) {
                    palette.base(back_base)
                } else {
                    palette.back()
                },
            )?;
            draw_base(
                out,
                y,
                front_x,
                front_base as char,
                if matches!(palette, Palette::Bases) {
                    palette.base(front_base)
                } else {
                    palette.front()
                },
            )?;
        }
        if speed > 0.0 {
            self.theta += 0.16 * speed;
            self.tick = self.tick.wrapping_add(speed.ceil() as usize);
        }
        status_line(out, height, mode, seed, speed, scale, focus, palette)
    }
}

#[cfg(feature = "codon-wheel")]
struct CodonState {
    theta: f32,
    tick: usize,
}

#[cfg(feature = "codon-wheel")]
impl CodonState {
    fn new(seed: u64) -> Self {
        Self {
            theta: (seed % 64) as f32 / 64.0 * TAU,
            tick: (seed % 64) as usize,
        }
    }

    fn draw(
        &mut self,
        out: &mut impl Write,
        seed: u64,
        speed: f32,
        scale: f32,
        focus: bool,
        palette: Palette,
    ) -> io::Result<()> {
        let (width, height) = terminal_size();
        let center_x = (width / 2).max(1) as i32;
        let center_y = (height / 2).max(1) as i32;
        let radius = (((width.min(height * 2) as f32) * 0.24) * scale).clamp(3.0, 28.0);
        let active = self.tick % 64;

        clear(out)?;
        write!(out, "\x1b[{center_y};{center_x}H\x1b[97m5'\x1b[0m")?;

        for idx in 0..64 {
            let first = idx / 16;
            let second = (idx / 4) % 4;
            let third = idx % 4;
            let codon = [RNA[first], RNA[second], RNA[third]];
            let angle = self.theta + (idx as f32 / 64.0) * TAU;
            let spoke_angle = self.theta + (first as f32 / 4.0) * TAU;
            let is_active = idx == active;

            draw_polar(
                out,
                center_x,
                center_y,
                radius * 0.38,
                spoke_angle,
                RNA[first] as char,
                if matches!(palette, Palette::Bases) && is_active {
                    palette.active_base(RNA[first])
                } else if matches!(palette, Palette::Bases) {
                    palette.base(RNA[first])
                } else {
                    palette.codon(0, is_active)
                },
            )?;
            draw_polar(
                out,
                center_x,
                center_y,
                radius * 0.56,
                angle,
                RNA[second] as char,
                if matches!(palette, Palette::Bases) && is_active {
                    palette.active_base(RNA[second])
                } else if matches!(palette, Palette::Bases) {
                    palette.base(RNA[second])
                } else {
                    palette.codon(1, is_active)
                },
            )?;
            draw_polar(
                out,
                center_x,
                center_y,
                radius,
                angle,
                RNA[third] as char,
                if matches!(palette, Palette::Bases) && is_active {
                    palette.active_base(RNA[third])
                } else if matches!(palette, Palette::Bases) {
                    palette.base(RNA[third])
                } else {
                    palette.codon(2, is_active)
                },
            )?;

            if is_active {
                draw_polar(
                    out,
                    center_x,
                    center_y,
                    radius + 5.0,
                    angle,
                    '*',
                    palette.accent(),
                )?;
                write!(
                    out,
                    "\x1b[{};2H\x1b[97m{}{}{} -> {}\x1b[0m",
                    height.saturating_sub(1),
                    codon[0] as char,
                    codon[1] as char,
                    codon[2] as char,
                    amino_acid(&codon)
                )?;
            }
        }

        if speed > 0.0 {
            self.theta += 0.045 * speed;
            self.tick = self.tick.wrapping_add(speed.ceil() as usize);
        }
        status_line(
            out,
            height,
            AppMode::Codon,
            seed,
            speed,
            scale,
            focus,
            palette,
        )
    }
}

fn status_line(
    out: &mut impl Write,
    height: u16,
    mode: AppMode,
    seed: u64,
    speed: f32,
    scale: f32,
    focus: bool,
    palette: Palette,
) -> io::Result<()> {
    if focus {
        return Ok(());
    }
    write!(
        out,
        "\x1b[{height};2H\x1b[90m{} | seed {} | speed {:.2}x | scale {:.2}x | color {} | left/right mode | up/down speed | +/- scale | c color | f focus | q exits\x1b[0m",
        mode.label(),
        seed,
        speed,
        scale,
        palette.label()
    )
}

fn alphabet_seed(alphabet: &[u8; 4]) -> u64 {
    if alphabet == RNA { 0xA0C6 } else { 0xA7C6 }
}

#[cfg(feature = "codon-wheel")]
fn draw_polar(
    out: &mut impl Write,
    center_x: i32,
    center_y: i32,
    radius: f32,
    angle: f32,
    ch: char,
    color: &str,
) -> io::Result<()> {
    let x = center_x + (angle.cos() * radius * 2.0) as i32;
    let y = center_y + (angle.sin() * radius) as i32;
    if x > 0 && y > 0 {
        write!(out, "\x1b[{y};{x}H{color}{ch}\x1b[0m")?;
    }
    Ok(())
}

#[cfg(feature = "codon-wheel")]
fn amino_acid(codon: &[u8; 3]) -> &'static str {
    match codon {
        b"UUU" | b"UUC" => "Phe (F)",
        b"UUA" | b"UUG" | b"CUU" | b"CUC" | b"CUA" | b"CUG" => "Leu (L)",
        b"UCU" | b"UCC" | b"UCA" | b"UCG" | b"AGU" | b"AGC" => "Ser (S)",
        b"UAU" | b"UAC" => "Tyr (Y)",
        b"UAA" | b"UAG" | b"UGA" => "Stop",
        b"UGU" | b"UGC" => "Cys (C)",
        b"UGG" => "Trp (W)",
        b"CCU" | b"CCC" | b"CCA" | b"CCG" => "Pro (P)",
        b"CAU" | b"CAC" => "His (H)",
        b"CAA" | b"CAG" => "Gln (Q)",
        b"CGU" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => "Arg (R)",
        b"AUU" | b"AUC" | b"AUA" => "Ile (I)",
        b"AUG" => "Met (M)",
        b"ACU" | b"ACC" | b"ACA" | b"ACG" => "Thr (T)",
        b"AAU" | b"AAC" => "Asn (N)",
        b"AAA" | b"AAG" => "Lys (K)",
        b"GUU" | b"GUC" | b"GUA" | b"GUG" => "Val (V)",
        b"GCU" | b"GCC" | b"GCA" | b"GCG" => "Ala (A)",
        b"GAU" | b"GAC" => "Asp (D)",
        b"GAA" | b"GAG" => "Glu (E)",
        b"GGU" | b"GGC" | b"GGA" | b"GGG" => "Gly (G)",
        _ => "?",
    }
}

fn draw_bridge(out: &mut impl Write, y: i32, a: i32, b: i32) -> io::Result<()> {
    let lo = a.min(b) + 1;
    let hi = a.max(b);
    for x in lo..hi {
        write!(out, "\x1b[{y};{x}H\x1b[90m-")?;
    }
    Ok(())
}

fn draw_base(out: &mut impl Write, y: i32, x: i32, base: char, color: &str) -> io::Result<()> {
    if x > 0 {
        write!(out, "\x1b[{y};{x}H{color}{base}\x1b[0m")?;
    }
    Ok(())
}

fn complement(base: u8, alphabet: &[u8; 4]) -> u8 {
    match base {
        b'A' if alphabet == RNA => b'U',
        b'A' => b'T',
        b'T' => b'A',
        b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base,
    }
}

fn coding_sequence(alphabet: &[u8; 4], codon_count: usize, rng: &mut Rng) -> Vec<u8> {
    let codons = if alphabet == RNA {
        RNA_CODONS
    } else {
        DNA_CODONS
    };
    let mut sequence = Vec::with_capacity(codon_count * 3);
    sequence.extend_from_slice(if alphabet == RNA { b"AUG" } else { b"ATG" });
    for _ in 1..codon_count {
        sequence.extend_from_slice(codons[rng.range(codons.len() as u32) as usize]);
    }
    sequence
}

fn clear(out: &mut impl Write) -> io::Result<()> {
    write!(out, "\x1b[H\x1b[2J")
}

fn sleep_frame(start: Instant, frame: Duration) {
    if let Some(remaining) = frame.checked_sub(start.elapsed()) {
        thread::sleep(remaining);
    }
}

fn scaled_frame(frame: Duration, speed: f32) -> Duration {
    if speed == 0.0 {
        return Duration::from_millis(80);
    }
    Duration::from_secs_f32(frame.as_secs_f32() / speed.max(0.01))
}

fn terminal_size() -> (u16, u16) {
    let output = Command::new("stty").arg("size").output();
    let Ok(output) = output else {
        return (80, 24);
    };
    let text = String::from_utf8_lossy(&output.stdout);
    let mut parts = text.split_whitespace();
    let height = parts.next().and_then(|s| s.parse().ok()).unwrap_or(24);
    let width = parts.next().and_then(|s| s.parse().ok()).unwrap_or(80);
    (width, height)
}

struct Column {
    head: i32,
    len: i32,
    speed: i32,
}

impl Column {
    fn reset(&mut self, height: i32, rng: &mut Rng) {
        self.head = -(rng.range(height.max(1) as u32) as i32);
        self.len = 4 + rng.range(18) as i32;
        self.speed = 1 + (rng.range(3) == 0) as i32;
    }
}

struct Columns(Vec<Column>);

impl Columns {
    fn new(width: usize, height: i32, rng: &mut Rng) -> Self {
        let mut cols = Vec::with_capacity(width);
        for _ in 0..width {
            let mut col = Column {
                head: 0,
                len: 0,
                speed: 1,
            };
            col.reset(height, rng);
            cols.push(col);
        }
        Self(cols)
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complements_dna() {
        assert_eq!(complement(b'A', DNA), b'T');
        assert_eq!(complement(b'T', DNA), b'A');
        assert_eq!(complement(b'C', DNA), b'G');
        assert_eq!(complement(b'G', DNA), b'C');
    }

    #[test]
    fn complements_rna() {
        assert_eq!(complement(b'A', RNA), b'U');
        assert_eq!(complement(b'U', RNA), b'A');
        assert_eq!(complement(b'C', RNA), b'G');
        assert_eq!(complement(b'G', RNA), b'C');
    }

    #[test]
    fn coding_sequence_keeps_frame_and_alphabet() {
        let mut rng = Rng::new(1);
        let dna = coding_sequence(DNA, 8, &mut rng);
        assert_eq!(&dna[..3], b"ATG");
        assert_eq!(dna.len(), 24);
        assert!(dna.iter().all(|base| DNA.contains(base)));

        let rna = coding_sequence(RNA, 8, &mut rng);
        assert_eq!(&rna[..3], b"AUG");
        assert_eq!(rna.len(), 24);
        assert!(rna.iter().all(|base| RNA.contains(base)));
    }

    #[test]
    fn cli_parses_short_aliases() {
        let args = ["m".to_string(), "r".to_string()];
        let cli = Cli::parse(&args).unwrap();
        assert_eq!(cli.mode, RequestedMode::Matrix);
        assert_eq!(cli.alphabet, Alphabet::Rna);

        let args = ["h".to_string(), "-d".to_string()];
        let cli = Cli::parse(&args).unwrap();
        assert_eq!(cli.mode, RequestedMode::Helix);
        assert_eq!(cli.alphabet, Alphabet::Dna);
    }

    #[test]
    fn cli_defaults_to_dna_helix() {
        let cli = Cli::parse(&[]).unwrap();
        assert_eq!(cli.mode, RequestedMode::Helix);
        assert_eq!(cli.alphabet, Alphabet::Dna);
    }

    #[test]
    fn cli_rejects_unknown_arguments() {
        let args = ["garbage".to_string()];
        assert!(Cli::parse(&args).is_err());
    }

    #[test]
    #[cfg(feature = "codon-wheel")]
    fn amino_acid_table_has_key_codons() {
        assert_eq!(amino_acid(b"AUG"), "Met (M)");
        assert_eq!(amino_acid(b"UGG"), "Trp (W)");
        assert_eq!(amino_acid(b"UAA"), "Stop");
        assert_eq!(amino_acid(b"GCU"), "Ala (A)");
    }
}
