# -*- coding: utf-8 -*-
"""Build the Scientific Reports submission manuscript (Chinese) from the
LetPub Scientific Reports template.

Content source:   manuscript_src/_norm.txt          (run-level normalised text)
Reference source: manuscript_src/_refs_nature.json   (formatted reference list)
Template:         manuscript_src/tpl.docx            (Scientific Reports Word template)
Figures:          figures/fig{1..6}_*.png
Output:           manuscript_out/_build.docx
"""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

import docx
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt

PROJ = Path(__file__).resolve().parent.parent.parent
SRC = PROJ / "manuscript_src"
TPL = SRC / "tpl.docx"
OUT_DIR = PROJ / "manuscript_out"
FIG_DIR = PROJ / "figures"
WORK = OUT_DIR / "_build.docx"
OUT_DIR.mkdir(parents=True, exist_ok=True)
raw = (SRC / "_norm.txt").read_text(encoding="utf-8")
blocks: dict[int, str] = {}
cur = None
for line in raw.split("\n"):
    m = re.match(r"### (\d+) \| (.+)", line)
    if m:
        cur = int(m.group(1))
        blocks[cur] = []
    elif cur is not None:
        blocks[cur].append(line)
P = {k: "\n".join(v).strip() for k, v in blocks.items()}

# ---- legacy figure call-outs -> the six-figure scheme -----------------------
FIG_REMAP = {
    "（图 2、图 3）": "（图 1、表 1）",
    "（图 4、图 5）": "（图 2）",
    "（图 7）": "（图 3）",
    "（图 9、图 10）": "（图 4）",
    "（图 11、图 12）": "（图 5）",
    "（图 13、图 14）": "（图 6a、b）",
    "（图 15 至图 17）": "（图 6c、表 2）",
}


def body(i: int) -> str:
    t = P[i]
    for a, b in FIG_REMAP.items():
        t = t.replace(a, b)
    # plain-text exponents in Methods -> superscript markers (keep the sign)
    t = re.sub(r"×10-(\d+)", r"×10^{-\1}", t)
    return t


TITLE = P[0]

ABSTRACT = P[2]
for lab in ("目的：", "方法：", "结果：", "结论："):
    ABSTRACT = ABSTRACT.replace(lab, "")
ABSTRACT = re.sub(r"\s+", "", ABSTRACT) if False else ABSTRACT
ABSTRACT = ABSTRACT.replace("本研究串联遗传因果、细胞状态与药理逆转。", "本研究串联遗传因果、细胞状态与药理逆转。")

KEYWORDS = P[4]

INTRO = [body(i) for i in (6, 7, 8, 9)]

RESULTS_LEAD = body(24)
RESULTS = [
    ("可成药基因组池的规模与药理构成", [body(26)]),
    ("全池 cis-eQTL 孟德尔随机化筛查", [body(28)]),
    ("共定位与多效性排除", [body(30)]),
    ("球结膜单细胞图谱中的状态轴与跨谱系重现", [body(32), body(33)]),
    ("连接图谱逆转筛选的统计效度", [body(35)]),
    ("全表型安全性判决", [body(37)]),
    ("五轴合成与候选优先级", [body(39)]),
]

DISCUSSION = [body(i) for i in range(41, 48)]

METHODS = [
    ("数据来源与可成药基因组池", [body(12)]),
    ("cis-eQTL 孟德尔随机化", [body(14)]),
    ("共定位与 SMR/HEIDI", [body(16)]),
    ("单细胞谱系内分层与跨谱系重现", [body(18)]),
    ("连接图谱逆转筛选", [body(20)]),
    ("全表型安全性扫描与五轴合成", [body(22)]),
]

# ------------------------------------------------------------------ references
refs = [r for r in json.loads((SRC / "_refs_nature.json").read_text(encoding="utf-8"))
        if re.match(r"^\[\d+\]", r)]
REF_FIX = {
    "[9] ": "[9] Chen, D. et al. Machine learning-driven identification of circulating "
            "biomarkers and therapeutic compounds for diabetic retinopathy. "
            "Hum. Genomics https://doi.org/10.1186/s40246-026-00967-2 (2026).",
}
REFS = []
for r in refs:
    n = int(re.match(r"^\[(\d+)\]", r).group(1))
    if n == 9:
        REFS.append(REF_FIX["[9] "])
    else:
        REFS.append(r)
assert len(REFS) == 32, len(REFS)

JOURNALS = [
    "Curr. Opin. Allergy Clin. Immunol.", "Front. Cell. Infect. Microbiol.",
    "Ann. Allergy Asthma Immunol.", "Transl. Vis. Sci. Technol.",
    "J. Allergy Clin. Immunol.", "Proc. Natl Acad. Sci. USA", "J. Immunol. Res.",
    "Cell Stem Cell", "Nat. Commun.", "Mol. Vis.", "Front. Immunol.",
    "Surv. Ophthalmol.", "Neuroimmunomodulation", "J. Headache Pain",
    "Atherosclerosis", "Acta Ophthalmol.", "Eye Vis.", "Int. J. Ophthalmol.",
    "PLoS Genet.", "Nat. Genet.", "Hum. Genomics", "Cell", "Nature", "Eye",
    "EBioMedicine", "Ocul. Surf.",
]
JOURNALS.sort(key=len, reverse=True)


# --------------------------------------------------------------- docx helpers
def clear_body(doc):
    b = doc.element.body
    for ch in list(b):
        if ch.tag == qn("w:sectPr"):
            continue
        b.remove(ch)


def add_par(doc, style, text=None, align=None, space_after=None, keep_next=False):
    p = doc.add_paragraph(style=style)
    if align is not None:
        p.alignment = align
    if space_after is not None:
        p.paragraph_format.space_after = space_after
    if keep_next:
        p.paragraph_format.keep_with_next = True
    if text:
        rich(p, text)
    return p


SUP_RE = re.compile(r"\^\{([^}]*)\}")


def rich(p, text, bold=False, italic=False):
    """Write text into paragraph, turning ^{...} into real superscript runs."""
    pos = 0
    for m in SUP_RE.finditer(text):
        if m.start() > pos:
            r = p.add_run(text[pos:m.start()])
            r.bold = bold or None
            r.italic = italic or None
        r = p.add_run(m.group(1))
        r.font.superscript = True
        r.bold = bold or None
        r.italic = italic or None
        pos = m.end()
    if pos < len(text):
        r = p.add_run(text[pos:])
        r.bold = bold or None
        r.italic = italic or None
    return p


def strip_style_numbering(doc, style_name):
    """The template's ReferenceList style auto-numbers each entry; the entries
    already carry their own [n] label, so drop the list numbering entirely."""
    el = doc.styles[style_name].element
    pPr = el.find(qn("w:pPr"))
    if pPr is not None:
        for old in pPr.findall(qn("w:numPr")):
            pPr.remove(old)


def ref_par(doc, text):
    """Reference entry with the journal name in italics and the volume in bold.
    Entries without a volume (advance-access articles carry only a DOI) still
    get the italic journal."""
    p = doc.add_paragraph(style="LetPub_7.2_ReferenceList")
    num, rest = re.match(r"^(\[\d+\]\s*)(.*)$", text, re.S).groups()
    p.add_run(num)

    hit = None
    for j in JOURNALS:
        idx = rest.find(j + " ")
        if idx != -1 and (hit is None or idx > hit[0]):
            hit = (idx, j)
    if hit is None:
        p.add_run(rest)
        return p
    idx, j = hit
    p.add_run(rest[:idx])
    r = p.add_run(j)
    r.italic = True
    after = rest[idx + len(j):]

    m = re.match(r"\s(\d+), ", after)
    if not m:
        p.add_run(after)
        return p
    p.add_run(" ")
    r = p.add_run(m.group(1))
    r.bold = True
    p.add_run(", " + after[m.end():])
    return p


def add_image(doc, path, width_mm=165):
    p = doc.add_paragraph(style="LetPub_4.1_Figure")
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(path), width=Mm(width_mm))
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(3)
    return p


def _border(el, edge, sz=8):
    tcPr = el.get_or_add_tcPr() if hasattr(el, "get_or_add_tcPr") else el
    borders = tcPr.find(qn("w:tcBorders"))
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tcPr.append(borders)
    e = borders.find(qn("w:" + edge))
    if e is None:
        e = OxmlElement("w:" + edge)
        borders.append(e)
    e.set(qn("w:val"), "single")
    e.set(qn("w:sz"), str(sz))
    e.set(qn("w:space"), "0")
    e.set(qn("w:color"), "000000")


def three_line(table):
    """Black-and-white three-line table: top rule, header rule, bottom rule."""
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    tblPr = table._tbl.tblPr
    for tag in ("w:tblBorders",):
        old = tblPr.find(qn(tag))
        if old is not None:
            tblPr.remove(old)
    b = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        e = OxmlElement("w:" + edge)
        e.set(qn("w:val"), "none")
        e.set(qn("w:sz"), "0")
        b.append(e)
    tblPr.append(b)

    for cell in table.rows[0].cells:
        _border(cell._tc, "top", 12)
        _border(cell._tc, "bottom", 6)
    for cell in table.rows[-1].cells:
        _border(cell._tc, "bottom", 12)


def fill_table(doc, header, rows, widths_mm, font_pt=9):
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.autofit = False
    tblPr = t._tbl.tblPr
    lay = OxmlElement("w:tblLayout")
    lay.set(qn("w:type"), "fixed")
    tblPr.append(lay)

    def cell_text(c, text, style_name, align):
        c.text = ""
        p = c.paragraphs[0]
        p.style = doc.styles[style_name]
        p.alignment = align
        pf = p.paragraph_format
        pf.line_spacing = 1.0
        pf.space_before = Pt(1)
        pf.space_after = Pt(1)
        r = p.add_run(str(text))
        r.font.size = Pt(font_pt)
        r.bold = style_name == "LetPub_5.2_TableColumnHeading"

    for j, h in enumerate(header):
        cell_text(t.rows[0].cells[j], h,
                  "LetPub_5.2_TableColumnHeading",
                  WD_ALIGN_PARAGRAPH.LEFT if j in (0, 1) else WD_ALIGN_PARAGRAPH.CENTER)
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            cell_text(t.rows[i].cells[j], v, "LetPub_5.2_TableBody",
                      WD_ALIGN_PARAGRAPH.LEFT if j in (0, 1) else WD_ALIGN_PARAGRAPH.CENTER)

    for row in t.rows:
        trPr = row._tr.get_or_add_trPr()
        cant = OxmlElement("w:cantSplit")
        trPr.append(cant)
        for j, w in enumerate(widths_mm):
            row.cells[j].width = Mm(w)
    three_line(t)
    return t


# --------------------------------------------------------------- figure texts
FIGURES = [
    ("fig1_study_design",
     "图 1　研究设计与分层可成药基因组池的构建。",
     "（a）六道关卡的分析流程：19,255 个蛋白编码基因与 DGIdb 的 26,739 个基因、5,011 个有药物交互记录的"
     "基因以及 19,935 个全血 cis-eQTL 数据集交叉，得到 10,983 个可成药基因，再经孟德尔随机化（6,329 个"
     "基因）、共定位与 SMR/HEIDI（47 个候选）、单细胞状态签名（600 个基因）、连接图谱逆转筛选（46 个"
     "化合物）与全表型安全性扫描（14 个基因 × 1,400 个端点）逐级收敛，最终给出 1 个 Tier 1 主候选。"
     "（b）三层药理证据的基因数与其中具备全血 cis-eQTL 汇总统计的基因数：T1 为 4,162 与 3,148、T2 为 "
     "5,028 与 3,621、T3 为 1,793 与 1,240。（c）各关卡保留的基因数，横轴为对数刻度，由 19,255 收敛至 1。"),
    ("fig2_mr_screen",
     "图 2　全池 cis-eQTL 孟德尔随机化筛查。",
     "（a）6,329 个基因 Wald ratio 的火山图，横轴为每 eQTL 等位基因的对数比值比，纵轴为负常用对数 P 值；"
     "三角形为 Bonferroni 显著者（2 个），实心圆为 FDR < 0.05 者（3 个），标注 IL18RAP、TLR1 与 TLR10。"
     "（b）六个基因的比值比与 95% 置信区间（实心圆），叠加独立样本中的重复估计（空心菱形）；粗体为三个"
     "先导基因：IL18RAP 0.92（0.89–0.94）、TLR1 1.42（1.24–1.63）、TLR10 1.57（1.28–1.93）。"
     "（c）一阶段 F 统计量的分布，中位 F = 347.7；629 个基因的 F 超过 3,000 未在图中显示，最大值为 40,299。"),
    ("fig3_coloc_smr",
     "图 3　共定位与 SMR/HEIDI 过滤。",
     "（a）50 个候选的 PP.H4 与 HEIDI 检验 P 值的对应关系，虚线为 PP.H4 = 0.8 与 HEIDI 检验 P = 0.05，"
     "标注 IL18RAP、TLR1、TLR10 与 NNT。（b）PP.H4 排序前 8 位候选的五种共定位后验概率分解，右侧数字为 "
     "PP.H4。（c）同一批候选的 SMR 显著性，虚线为 0.05/8 的 Bonferroni 阈值，右侧标注 HEIDI 检验是否通过；"
     "TLR1 的 SMR P = 9.38×10^{-7}、HEIDI P = 0.092。"),
    ("fig4_singlecell",
     "图 4　球结膜单细胞图谱中的表达状态轴与跨谱系重现。",
     "（a）TLR1 在八个球结膜谱系中的检出率，范围 11.6% 至 20.6%；深色为进入主分析的谱系，浅色为每组"
     "细胞数不足 200 而被排除的上皮谱系，柱内为检出细胞数与总细胞数。（b）按靶基因表达切分的两组细胞"
     "在总 UMI、TLR1 表达量与残差上的平衡，总 UMI 差 −3.00%，TLR1 差 +57.68 CPM，残差差 1.02 s.d.。"
     "（c）600 个签名基因的对数倍数变化与显著性，标注上调臂 ITM2A、PDLIM7、SPARCL1、BACH2 与下调臂 "
     "SAA1、CD200R1。（d）签名在六个谱系中的重现，实线为与主分析的 Pearson 相关系数，点线为全基因符号"
     "一致率，两条水平虚线为各自的随机分层阴性对照；平均方向一致率 73.5%，阴性对照 r = +0.0101。"),
    ("fig5_cmap",
     "图 5　连接图谱逆转筛选的统计效度。",
     "（a）四个查询规模下观测到的最高连接分与同规模随机 landmark 集零分布的对照（柱为均值，误差棒为 "
     "1 s.d.），标注 z 值与置换 P 值；15 个基因时 z = 2.47、P = 0.0476。（b）逆转特异性排序前 16 位的"
     "化合物，BI-2536 为 0.0949，其余 15 个并列 0.0324。（c）糖皮质激素阳性对照在四个查询规模下的最佳"
     "排名，13 个成员中最高排名为 1 位。"),
    ("fig6_safety_integration",
     "图 6　全表型安全性判决与五轴合成。",
     "（a）14 个候选基因在 1,400 个 FinnGen 端点中的显著关联数，按疾病域分层；TLR1 为 12 个、TLR10 为 "
     "14 个、PBXIP1 为 66 个，右侧数字为眼科信号数与重症信号数。（b）TLR1 与 6 个 Tier 2 支持性候选在"
     "五个轴上的得分，空白格表示该基因未进入对应分析（得分为 0）。（c）综合分前 10 位候选，虚线为 0.55 "
     "的 Tier 1 阈值，上方括号标注候选层级。"),
]

# --------------------------------------------------------------- declarations
ETHICS = ("本研究仅使用已公开发表的汇总统计量与公开单细胞数据集，不涉及新的生物样本采集、"
          "受试者招募或动物实验，无需伦理委员会审批。所有原始数据的采集均遵循原研究已获批准的"
          "伦理规程与《赫尔辛基宣言》。")

DATA_AVAIL = ("本研究分析的全部数据均来自公开资源。单细胞转录组数据来自 Gene Expression Omnibus "
              "（GSE242382）。过敏性结膜炎的全基因组关联汇总统计来自 FinnGen 第 13 版"
              "（https://www.finngen.fi/）。全血 cis-eQTL 汇总统计来自 eQTLGen 联盟，并通过 IEU OpenGWAS "
              "项目获取（https://gwas.mrcieu.ac.uk/）。药物与靶点注释来自 DGIdb 与 HGNC，扰动参考谱来自 "
              "LINCS L1000。本研究未产生新的原始数据；全部中间结果表与分析日志已随本文一并提供。")

CODE_AVAIL = ("用于可成药基因组池构建、孟德尔随机化筛查、共定位与 SMR/HEIDI 分析、单细胞谱系内分层、"
              "连接图谱逆转查询与全表型扫描的自编代码托管于 https://github.com/[USER]/[REPO]，"
              "并已归档至 Zenodo（https://doi.org/10.5281/zenodo.[XXXXXXX]）。分析阈值固化于 "
              "config/config.json。其余分析使用公开软件包完成，版本见方法部分。")

ACK = ("感谢 FinnGen 研究团队、eQTLGen 联盟与 IEU OpenGWAS 项目公开共享汇总统计资源，"
       "感谢 GSE242382 数据集的作者公开球结膜单细胞转录组数据。")

AUTHOR_CONTRIB = ("[AUTHOR 1] 构思并设计研究，完成数据整理、质量控制与下游分析，撰写论文初稿。"
                  "[AUTHOR 2] 完成孟德尔随机化分析、共定位与敏感性检验。"
                  "[AUTHOR 3] 完成单细胞分析与全部图件制作。"
                  "全体作者参与结果解读与论文修改，并批准最终版本。")

FUNDING = "本研究由 [FUNDING AGENCY] 资助，项目编号 [GRANT NUMBER]。"

COMPETING = "全体作者声明无利益冲突。"

STROBE = ("本文按 STROBE-MR 清单逐条报告：工具变量与暴露、结局的对应关系，样本重叠的界定，"
          "弱工具偏倚的检验，以及水平多效性的多重判据，均在方法各小节中给出对应条目与阈值。"
          "完整清单随投稿一并提交。")

AI_DECL = ("生成式人工智能工具被用于三项辅助性工作：分析代码的初稿编写、图件版式与配色的组装、"
           "以及中文表述的语言润色。所有由人工智能辅助产出的代码均经作者逐行审阅并在完整数据上"
           "重新运行验证，全部图件均由作者依据本文报告的数值复核，文字内容由作者改写与定稿。"
           "人工智能工具未参与研究设计、结果解释或结论形成，作者对本文的全部内容承担责任。")

# ---------------------------------------------------------------------- build
OUT_DIR.mkdir(parents=True, exist_ok=True)
shutil.copy(TPL, WORK)
doc = docx.Document(str(WORK))
clear_body(doc)
strip_style_numbering(doc, "LetPub_7.2_ReferenceList")
for _st in ("LetPub_4.2_FigureCaption", "LetPub_5.1_TableCaption",
            "LetPub_5.3_TableFootnote", "LetPub_4.1_Figure"):
    _pf = doc.styles[_st].paragraph_format
    _pf.line_spacing = 1.0
    _pf.space_before = Pt(4)
    _pf.space_after = Pt(2)

# title block
add_par(doc, "LetPub_1.1_ArticleTitle", TITLE, align=WD_ALIGN_PARAGRAPH.CENTER)

p = doc.add_paragraph(style="LetPub_1.2_Author")
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.add_run("[AUTHOR 1]")
r = p.add_run("1,+"); r.font.superscript = True
p.add_run(", [AUTHOR 2]")
r = p.add_run("2,+"); r.font.superscript = True
p.add_run(" & [AUTHOR 3]")
r = p.add_run("1"); r.font.superscript = True

for n in ("1", "2", "3"):
    q = doc.add_paragraph(style="LetPub_1.3_Affiliation")
    r = q.add_run(n); r.font.superscript = True
    q.add_run("[DEPARTMENT, INSTITUTION, CITY, POSTCODE, COUNTRY]")

q = doc.add_paragraph(style="LetPub_1.4_Correspondence")
q.add_run("*e-mail: [CORRESPONDING AUTHOR E-MAIL]")
q = doc.add_paragraph(style="LetPub_1.4_Correspondence")
q.add_run("+共同第一作者。")

# abstract
add_par(doc, "LetPub_1.5_AbstractHead", "摘要")
add_par(doc, "LetPub_1.6_AbstractText", ABSTRACT)
add_par(doc, "LetPub_1.7_KeywordsHead", "关键词")
add_par(doc, "LetPub_1.8_Keywords", KEYWORDS)

# introduction
add_par(doc, "LetPub_2.1_Heading1", "引言")
for t in INTRO:
    add_par(doc, "LetPub_3.2_TextIndent", t)

# results
add_par(doc, "LetPub_2.1_Heading1", "结果")
add_par(doc, "LetPub_3.2_TextIndent", RESULTS_LEAD)
for head, paras in RESULTS:
    add_par(doc, "LetPub_2.2_Heading2", head)
    for t in paras:
        add_par(doc, "LetPub_3.2_TextIndent", t)

# discussion
add_par(doc, "LetPub_2.1_Heading1", "讨论")
for t in DISCUSSION:
    add_par(doc, "LetPub_3.2_TextIndent", t)

# methods
add_par(doc, "LetPub_2.1_Heading1", "方法")
for head, paras in METHODS:
    add_par(doc, "LetPub_2.2_Heading2", head)
    for t in paras:
        add_par(doc, "LetPub_3.2_TextIndent", t)
add_par(doc, "LetPub_2.3_Heading3", "报告规范")
add_par(doc, "LetPub_3.2_TextIndent", STROBE)

add_par(doc, "LetPub_2.3_Heading3", "人工智能使用声明")
add_par(doc, "LetPub_3.2_TextIndent", AI_DECL)

# ethics
add_par(doc, "LetPub_2.2_Heading2", "伦理声明")
add_par(doc, "LetPub_3.2_TextIndent", ETHICS)

# references
add_par(doc, "LetPub_7.1_References", "参考文献")
for r_ in REFS:
    ref_par(doc, r_)

# back matter, in the order Scientific Reports specifies
add_par(doc, "LetPub_6.1_Declarations", "致谢")
add_par(doc, "LetPub_3.1_Text", ACK)
add_par(doc, "LetPub_6.1_Declarations", "作者贡献")
add_par(doc, "LetPub_3.1_Text", AUTHOR_CONTRIB)
add_par(doc, "LetPub_6.1_Declarations", "基金资助")
add_par(doc, "LetPub_3.1_Text", FUNDING)
add_par(doc, "LetPub_6.1_Declarations", "数据可用性")
add_par(doc, "LetPub_3.1_Text", DATA_AVAIL)
add_par(doc, "LetPub_6.1_Declarations", "代码可用性")
add_par(doc, "LetPub_3.1_Text", CODE_AVAIL)
add_par(doc, "LetPub_6.1_Declarations", "利益冲突")
add_par(doc, "LetPub_3.1_Text", COMPETING)

# figure legends
add_par(doc, "LetPub_2.1_Heading1", "图注")
for key, title, desc in FIGURES:
    p = doc.add_paragraph(style="LetPub_4.2_FigureCaption")
    rich(p, title, bold=True)
    rich(p, desc)

# tables
add_par(doc, "LetPub_2.1_Heading1", "表格")

add_par(doc, "LetPub_5.1_TableCaption",
        "表 1　分层可成药基因组池的构成与 cis-eQTL 工具化情况。", keep_next=True)
fill_table(
    doc,
    ["层级", "定义", "基因数", "有全血 cis-eQTL", "占比"],
    [
        ["T1", "有临床期注释且有药物交互", "4,162", "3,148", "75.6%"],
        ["T2", "仅有临床期注释", "5,028", "3,621", "72.0%"],
        ["T3", "仅有文献或实验证据", "1,793", "1,240", "69.2%"],
        ["合计", "分层可成药池", "10,983", "8,009", "72.9%"],
    ],
    [16, 58, 22, 46, 23],
)
add_par(doc, "LetPub_5.3_TableFootnote",
        "占比为各层具备全血 cis-eQTL 汇总统计的基因占该层基因数的百分比。T1 与 T2 合计占全池的 83.6%。")

doc.add_page_break()
add_par(doc, "LetPub_5.1_TableCaption",
        "表 2　五轴合成排序前 7 位候选基因的得分构成。", keep_next=True)
fill_table(
    doc,
    ["排名", "基因", "MR", "共定位", "单细胞", "CMap", "安全性", "综合分", "候选层级"],
    [
        ["1", "TLR1", "0.8270", "0.9776", "0.3733", "0.9493", "0.55", "0.7631", "Tier 1"],
        ["2", "TXNDC12", "0.4269", "0.7537", "0", "0", "1.00", "0.3788", "Tier 2"],
        ["3", "PRKAB1", "0.3131", "0.9203", "0", "0", "1.00", "0.3780", "Tier 2"],
        ["4", "TLR10", "0.5685", "0.7235", "0", "0", "0.55", "0.3703", "Tier 2"],
        ["5", "PHF1", "0.5036", "0.5925", "0", "0", "1.00", "0.3696", "Tier 2"],
        ["6", "IL18RAP", "1.0000", "0.1337", "0", "0", "0.40", "0.3667", "Tier 2"],
        ["7", "ATP7B", "0.3005", "0.8136", "0", "0", "1.00", "0.3529", "Tier 2"],
    ],
    [15, 24, 16, 17, 17, 15, 17, 18, 26],
)
add_par(doc, "LetPub_5.3_TableFootnote",
        "综合分按固定权重合成：MR 0.3、共定位 0.2、单细胞 0.2、CMap 0.2、安全性 0.1，"
        "0.55 以上判为 Tier 1。单细胞与 CMap 得分为 0 表示该基因未进入对应层级的分析。")

# figures
doc.add_page_break()
add_par(doc, "LetPub_2.1_Heading1", "图件")
missing = []
for key, title, desc in FIGURES:
    f = FIG_DIR / (key + ".png")
    if not f.exists():
        missing.append(str(f))
        continue
    add_image(doc, f)
    p = doc.add_paragraph(style="LetPub_4.2_FigureCaption")
    rich(p, title, bold=True)

doc.save(str(WORK))
print("saved", WORK)
if missing:
    print("MISSING FIGURES:", *missing, sep="\n  ")