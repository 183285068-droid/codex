# Codex 当前云端基础运行环境

已在本次云端会话配置并验证 Python 数据处理、基础统计、机器学习和绘图环境。
这是供当前助手执行代码的环境。没有额外安装独立 Codex CLI，也没有改动宿主模型、账号、审批或网络策略。

## 使用

在本目录运行：

```bash
bash setup.sh
./run.sh scripts/check_env.py
./run.sh your_script.py
```

先进入克隆仓库中的 `codex-environment/` 目录，再执行上述命令。
运行脚本会先进入该目录，文件参数请相对该目录填写，或使用绝对路径。
数据可放在 `data/`，结果可写入 `outputs/`。不建议把唯一一份原始数据只存放在临时云端目录。

## 已配置内容

- Python 3.12.14；项目虚拟环境 `.venv` 复用宿主已安装的软件包，新增安装限制在虚拟环境内。
- NumPy、pandas、SciPy、Matplotlib、seaborn、scikit-learn、openpyxl，具体版本见 requirements.txt。
- UTF-8 文本处理、无界面绘图、项目内缓存目录。
- 自检覆盖 CSV/Excel 读写、Welch t 检验、已知答案的线性回归、PNG 图形生成。
- `outputs/environment-check/` 是自检使用的人工示例数据及结果，并非实际科研分析结果。
- 宿主 Node.js 24.19.0、npm 11.9.0、Git 2.51.1 已可用；未另建 JavaScript 项目。

## 恢复与限制

云端运行目录可能被回收；解压此配置包后执行 `bash setup.sh` 可重新配置。
若宿主已包含匹配版本的软件包，设置不需要联网；否则会尝试从当前 Python 软件源安装。
安装能否成功取决于宿主网络权限和软件源可用性。requirements.txt 固定直接依赖，不是跨平台的完整传递依赖锁文件。
配置包不包含虚拟环境、软件包本体、用户数据、密钥或登录信息。

当前未安装 R、Bioconductor、lifelines、statsmodels 或 Jupyter。
因此此配置并不表示已具备 DESeq2、GO/KEGG 富集或 Cox 回归的完整专用环境；需要具体分析时再添加相应工具。
也不代表后续所有新会话都会自动采用此目录；新会话可提供配置包并要求恢复。

## 官方配置说明

Codex 客户端配置和这里的 Python 代码执行环境属于不同配置层。
官方说明：https://learn.chatgpt.com/docs/config-file/config-basic
