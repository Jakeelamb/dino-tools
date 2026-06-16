//! Hand-written x86-64 AVX2 newline mask helper for [`super::scan_newlines`].
//!
//! Enabled with the `asm-scan` crate feature. Intended for experiments comparing
//! LLVM-generated intrinsics against a fixed assembly sequence.

use std::arch::global_asm;

global_asm!(
    ".text",
    ".globl dino_seq_avx2_newline_mask_32",
    ".p2align 4",
    "dino_seq_avx2_newline_mask_32:",
    "    vpbroadcastb ymm0, byte ptr [rip + newline_byte]",
    "    vmovdqu ymm1, [rdi]",
    "    vpcmpeqb ymm1, ymm1, ymm0",
    "    vpmovmskb eax, ymm1",
    "    vzeroupper",
    "    ret",
    ".section .rodata",
    ".p2align 1",
    "newline_byte:",
    "    .byte 10",
    ".section .text",
);

unsafe extern "C" {
    fn dino_seq_avx2_newline_mask_32(ptr: *const u8) -> u32;
}

#[target_feature(enable = "avx2")]
unsafe fn scan_newlines_avx2_asm_inner(bytes: &[u8], out: &mut Vec<usize>) {
    const LANES: usize = 32;
    let ptr = bytes.as_ptr();
    let mut i = 0usize;
    while i + LANES <= bytes.len() {
        let mask = unsafe { dino_seq_avx2_newline_mask_32(ptr.add(i)) };
        let mut mask = mask;
        while mask != 0 {
            let bit = mask.trailing_zeros() as usize;
            out.push(i + bit);
            mask &= mask - 1;
        }
        i += LANES;
    }
    out.extend(memchr::memchr_iter(b'\n', &bytes[i..]).map(|offset| i + offset));
}

pub(super) fn scan_newlines_impl(bytes: &[u8], out: &mut Vec<usize>) {
    if std::is_x86_feature_detected!("avx2") {
        unsafe {
            scan_newlines_avx2_asm_inner(bytes, out);
        }
    } else {
        out.extend(memchr::memchr_iter(b'\n', bytes));
    }
}
