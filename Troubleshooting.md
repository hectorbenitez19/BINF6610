I used ChatGpt assistance when writing this assignment. I used it when I would get stuck and troubleshoot what the issue is and what steps I need to take to resolve it.

I ran into an issue where I was failing on my check for if the reference dictionary exists or not in stage 0. I was stumped because my previous test to check for the fa file in the reference fai file exists and it does. I found that the issue was that I never added it as a local variable. Since the ref fa file was already set as a global variable my fai test would pass and check properly so I added a variable for the dict and it passed.


In stage 4 I ran into an issue when running the stage [E::hts_open_format] Failed to open file "/home/hector/smoke-out/postprocess/smoke_01.sorted.bam" : No such file or directory
samtools sort: failed to create "/home/hector/smoke-out/postprocess/smoke_01.sorted.bam": No such file or directory
The issue is that the directory where the sorted bam files would be stored didnt get created yet before the stage tries to access it. The cause of this problem is due to me not adding my POST="${OUT}/postprocess" with the other directories so it would not get added. This was a syntax issue and I was able to resolve it by following the structure of how the other directories were getting added together.