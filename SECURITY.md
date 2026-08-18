# Security Policy

## Supported versions

Security fixes target the current `main` branch and the latest accepted release
named by `governance/releases/current.yaml`. Older snapshots are immutable
provenance and may not receive backports.

## Reporting a vulnerability

Please use [GitHub private vulnerability
reporting](https://github.com/vantasnerdan/substrate-framework/security/advisories/new).
Do not open a public issue, discussion, or pull request for an undisclosed
vulnerability, and do not include live credentials or unnecessary personal
data in a report.

Include:

- the affected commit, file, symbol, or workflow;
- the attack prerequisites and realistic impact;
- a minimal reproduction or proof of concept;
- whether the issue is already public or actively exploited; and
- any proposed mitigation, if known.

Maintainers will acknowledge a report, assess scope and severity, coordinate a
fix and disclosure plan, and credit reporters who want attribution. Please
allow a reasonable remediation period before public disclosure.

## Security scope

Security reports include credential exposure, dependency vulnerabilities,
unsafe deserialization or code execution, workflow privilege escalation,
untrusted pull-request secret access, and integrity defects in release or
governance tooling. A disputed scientific conclusion without a software or
process-security impact belongs in the normal issue and claim-review workflow.

Never use a vulnerability report to submit confidential third-party research
material that the repository lacks rights to redistribute.
