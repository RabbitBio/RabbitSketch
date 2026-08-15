#!/usr/bin/env python3
"""Independent Python check for the RankStream v1 C++ golden vectors."""

MASK64 = (1 << 64) - 1
SEED = 42

ADMISSION = 0x243F6A8885A308D3
STEP = 0x9E3779B97F4A7C15
DERIVE_SEED = 0xD1B54A32D192ED03

VECTORS = (
    (0, 0x810879608E4259CC, 0x178C233B2E1AEC71),
    (1, 0x203EA4C5049AD615, 0x77F9F74BF0D207EC),
    (6, 0x3B7577DA105E355B, 0x93D5850733F0AD06),
    (27, 0xCF8FFB89367B9DB1, 0xB7CCB3C7C6587849),
    (
        116_418_878_235,
        0xCDBA8A2A8AA4D972,
        0x358F8FA87D4B46DF,
    ),
)


def fmix64(value: int, seed: int) -> int:
    value = (value ^ seed) & MASK64
    value ^= value >> 33
    value = (value * 0xFF51AFD7ED558CCD) & MASK64
    value ^= value >> 33
    value = (value * 0xC4CEB9FE1A85EC53) & MASK64
    value ^= value >> 33
    return value & MASK64


def derive(fingerprint: int, domain: int, counter: int = 0) -> int:
    step = (STEP * (counter + 1)) & MASK64
    return fmix64(
        fingerprint ^ domain ^ step,
        SEED ^ DERIVE_SEED,
    )


def main() -> int:
    for canonical, expected_fp, expected_admission in VECTORS:
        fingerprint = fmix64(canonical, SEED)
        assert fingerprint == expected_fp
        assert derive(fingerprint, ADMISSION) == expected_admission
    print("test_rankstream_vectors.py: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
