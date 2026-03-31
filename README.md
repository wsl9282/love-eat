# love-eat

## 系统发育树脚本（Python）

新增脚本：`phylogenetic_tree.py`，用于从**已比对（等长）**的 FASTA 序列构建 UPGMA 系统发育树，并输出 Newick 格式。

### 用法

```bash
python3 phylogenetic_tree.py -i input.fasta -o tree.nwk
```

不指定 `-o` 时会直接打印到标准输出：

```bash
python3 phylogenetic_tree.py -i input.fasta
```
