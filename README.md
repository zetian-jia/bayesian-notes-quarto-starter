# 贝叶斯数据分析笔记

跟着 Q1–Q5 框架整理吴喜之《贝叶斯数据分析——基于R与Python的实现》的阅读笔记，用 [Quarto](https://quarto.org) 构建成静态网站。

## 本地预览

```bash
# 安装 Quarto: https://quarto.org/docs/get-started/
quarto preview
```

## 目录结构

```
bayesian-notes/
├── _quarto.yml          # book 项目配置：章节顺序、主题
├── index.qmd            # 首页/前言
├── styles.css           # 自定义样式
├── chapters/
│   ├── 01-intro.qmd
│   ├── 02-basic-concepts.qmd
│   ├── 03-bernoulli.qmd
│   ├── 04-poisson.qmd
│   ├── 05-normal.qmd
│   ├── 06-mcmc-algorithms.qmd
│   └── 07-stan-pymc.qmd
└── .github/workflows/publish.yml   # push 到 main 自动构建部署
```

## 部署

push 到 `main` 分支会自动触发 GitHub Actions，渲染后发布到 `gh-pages` 分支。
第一次使用前，去仓库 **Settings → Pages**，把 Source 设为 `gh-pages` 分支。

## 新增一章

1. 在 `chapters/` 下新建 `.qmd` 文件
2. 在 `_quarto.yml` 的 `book.chapters` 里加一行引用它
3. push，自动重新构建
