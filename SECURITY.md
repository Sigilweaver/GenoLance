# Security policy

## Supported versions

Only the latest minor release of GenoLance is supported with
security fixes.

| Version | Supported |
| ------- | --------- |
| 0.2.x   | Yes       |
| < 0.2   | No        |

## Reporting a vulnerability

Please report security vulnerabilities privately via
[GitHub Security Advisories](https://github.com/Sigilweaver/GenoLance/security/advisories/new).

Do **not** open a public issue for security reports. We will
acknowledge within 7 days and aim to publish a fix or mitigation
within 30 days for confirmed issues.

## Scope

In scope:

- Memory-safety bugs in the variant-store parser or ingest path.
- Decompression bombs, path traversal, or arbitrary file write
  triggered by a malformed `.vcf`, `.vcf.gz`, or Lance dataset.
- Crashes triggered by an attacker-supplied file.
- Supply-chain integrity issues affecting published crates.

Out of scope:

- Incorrect variant interpretation. GenoLance is not validated for
  clinical use; open a normal issue with a reproducer.
- Vulnerabilities in upstream Lance / Arrow / htslib crates:
  please report those upstream.

## Disclosure

Coordinated disclosure is preferred. Once a fix is released, the
advisory will be made public and credited to the reporter unless
they request anonymity.
