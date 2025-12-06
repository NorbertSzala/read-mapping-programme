# read-mapping-programme

# Algorithms for Genomic Data Analysis  
## Assignment 1: Read Mapping  
**Winter semester 2025/2026**

---

## **Task**

Implement a read mapping algorithm that:

- is designed to work on reads of length ∼1 kbp with error rate **5–10%**,
- may fail to map some reads, but **should avoid incorrect alignments**,
- is **efficient** (fast) and of **good quality** (high proportion of mapped reads).

Your program can be based on any approach to the **approximate string matching** problem.

In your program you may use:

- code fragments from classes,
- libraries included in standard Python distribution,
- **NumPy, SciPy, Biopython**,
- code for building suffix arrays (e.g., **Karkkainen–Sanders algorithm**).

You **cannot** use:

- programs/libraries for assembly, mapping, alignment, etc.,
- multiprocessing commands,
- subprograms written in other languages,
- JIT compilers.

The solution should include:

- a **Python 3** program file,
- slides with a short description of your approach.

Avoid submitting a packaged virtual environment.

---

## **Specification**

Minimum performance requirements:

Given a reference sequence of length **≤ 20 Mbp**, a collection of *r* reads must be processed in:

- **≤ 5 + r/10 minutes** (students' server),
- **< 1 GB RAM**,
- **≥ 80%** of mappable reads correctly mapped,
- **≤ 1%** mapped incorrectly.

Program must be executable as:

python3 mapper.py reference.fasta reads.fasta output.txt

**Input**: FASTA files  
**Output**: For each mapped read, one line:
<read_id> <start> <end>

(tab-separated)

---

## **Assumptions**

- Most reads come from the reference but contain errors (SNPs, insertions, deletions).  
- Errors occur independently at 5–10%, total number may slightly exceed 10%.  
- Input may include **unmappable reads** (error rate > 20%).  
- A read is **correctly mapped** if:
  - it is mappable, and  
  - reported coordinates differ by **≤ 20 bp** from the true source.

---

## **Attached Files**

**Package 1:**
- Python implementation of Karkkainen–Sanders SA algorithm,
- example mapping program (not sufficient for minimal requirements),
- small example input/output files.

**Package 2:**
- full-size example input/output files.

---

## **Terms and Conditions**

Can be completed **individually** or in **2–3 person teams**.

### **Schedule**
- Team submission to: ******* — **November 4**
- Solution submission to Moodle — **November 23**
- In-class presentation — **December 16**

---

## **Assessment**

Solutions meeting minimal requirements → **2 points**.

Additional points:

### **Mapping correctness**
- 3 pts — **0%** incorrectly mapped reads  
- 2 pts — **≤ 0.2%** incorrectly mapped  
- 1 pt — **≤ 0.5%** incorrectly mapped  

### **Mapping completeness**
- 3 pts — **≥ 99%** reads mapped  
- 2 pts — **≥ 95%**  
- 1 pt — **≥ 90%**  

### **Mapping time**
- 3 pts — **≤ 1 sec/read**  
- 2 pts — **≤ 2 sec/read**  
- 1 pt — **≤ 3 sec/read**  

### **Other**
- up to **2 pts** — deadlines & presentation quality  
- **Team size:**  
  - 2 pts — 1 person  
  - 1 pt — 2 persons

---
