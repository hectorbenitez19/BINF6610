I used AI to help assist me when writing this assignment. The AI was used to help me troubleshoot when I get stuck on problems and refactor certian parts of the code.

1. During the laptop smoke test, Snakemake was executed with --forcerun haplotypecaller. The haplotypecaller job for smoke_02 completed successfully. While haplotypecaller was processing smoke_03, the workflow was interrupted using Ctrl+C. Snakemake displayed Terminating processes on user request, followed by a CalledProcessError and a WorkflowError indicating that at least one job did not complete successfully. The interruption occurred during an active job, demonstrating that Snakemake responds to user-requested termination.

The next step was to restart the workflow and verify that interrupted work could be recovered without unnecessarily rerunning completed upstream stages.

After interrupting the haplotypecaller job for smoke_03, I performed a dry run without --forcerun. Snakemake detected that work/gvcf/smoke_03.g.vcf.gz and its index were missing and scheduled the job for re-execution.

The completed smoke_01 and smoke_02 HaplotypeCaller jobs were not scheduled again. Earlier alignment and duplicate-marking stages were also preserved. Snakemake scheduled the downstream joint-genotyping, filtering, MultiQC, and publishing steps as needed.

This demonstrates that Snakemake can identify interrupted work and plan recovery without restarting the entire pipeline. 

2. I investigated the failed haplotypecaller job from the interrupted laptop smoke-test run.

Diagnostic command: tail -n 45 .snakemake/log/2026-10-16T202313.734789.snakemake.log

The workflow log identified haplotypecaller as the failed rule, with job ID 12 and sample wildcard smoke_03. It displayed the GATK command, expected output files, and job-specific log path (logs/haplotypecaller/smoke_03.log). The error was reported as CalledProcessError with a nonzero exit status. The message Terminating processes on user request immediately preceded the failure, confirming that the error was associated with the intentional Ctrl+C interruption. I restarted the workflow without --forcerun. Snakemake regenerated the missing smoke_03 outputs and successfully completed the remaining pipeline stages.

 When a Snakemake job fails, inspect the workflow log in .snakemake/log/ to identify the failed rule, sample, command, and error. Then inspect the rule-specific log under logs/ for application-level details. Rule-specific logs may be overwritten by subsequent executions, so preserve them before retrying when investigating an unexpected failure.

 3. I Requested a FastQC output using an invalid read number (R3) to verify that the workflow rejects unsupported wildcard values.

Diagnostic command:

snakemake --profile profiles/laptop -n \
    work/fastqc/smoke_01_R3_fastqc.zip \
    --config samplesheet=$HOME/BINF6610/Assignment1/smoke/samplesheet.csv \
             ref=$HOME/BINF6610/Assignment1/smoke/smoke.fa \
             region=smoke_1mb

Snakemake raised MissingRuleException: No rule to produce work/fastqc/smoke_01_R3_fastqc.zip during DAG construction. The fastqc rule defines wildcard_constraints: read="[12]", restricting the read wildcard to 1 or 2. The requested R3 output did not satisfy this constraint. Use a valid read number (R1 or R2) and verify that the requested sample and read combination exists in the samplesheet. Wildcard constraints prevent invalid filenames from matching workflow rules. A MissingRuleException can indicate an invalid requested target or a wildcard mismatch, rather than a failure of the underlying bioinformatics tool.

4. I Checked whether Apptainer was available in the local WSL Ubuntu environment used for laptop smoke testing.

Diagnostic command:

command -v apptainer || echo "Apptainer is not installed or not on PATH"

Apptainer is not installed or not on PATH. Apptainer is not available through the current local environment's PATH. The Explorer profile specifies software-deployment-method: apptainer, so that profile requires Apptainer to be available in its execution environment. Use profiles/laptop for local smoke testing, where the workflow runs using locally installed bioinformatics tools. Use profiles/explorer on the Explorer HPC cluster, where the necessary Slurm and Apptainer environment must be configured and verified. A workflow may run successfully with one execution profile but fail with another because the environments have different software dependencies. Before running on a new platform, verify that the required container runtime, scheduler, and software environment are available.