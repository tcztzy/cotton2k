# Cotton2K GUI（Tauri v2）Proposal

## 1. 目标与范围

### 目标
- 为 `cotton2k` 提供跨平台桌面 GUI（macOS / Windows / Linux）。
- 保持现有 Rust 模拟内核为唯一计算真源，不重复实现算法。
- 优先交付“可跑、可看、可导出”的科研工作流：配置 -> 运行 -> 可视化 -> 归档。

### MVP 范围（第一阶段）
- 选择并加载 `profile.toml` / `profile.json`。
- 一键运行模拟（调用 Rust 内核）。
- 展示运行状态与错误信息。
- 读取并展示 `output.csv` 关键曲线（如 `lint_yield`、`leaf_area_index`、`plant_height`）。
- 导出结果目录（含输入快照 + 输出 CSV + 运行元数据）。

### 非目标（MVP 暂不做）
- 复杂 3D 可视化。
- 多任务分布式调度。
- 移动端（Tauri v2 支持移动端，但本提案聚焦桌面端）。

## 2. 调研结论（针对 Tauri v2）

- Tauri v2 在 2024-10-02 发布 stable，且后续已持续迭代（当前按 v2 最新小版本跟进），架构可用于生产。
- Tauri 的前后端通信模型适合本项目：前端 Web UI + Rust 命令桥接（`#[tauri::command]`）。
- v2 的权限模型更细：通过 `permissions` + `capabilities` 显式授权，可把文件访问、子进程执行限制到最小范围。
- 若需要隔离执行，可用 sidecar 机制把已有 CLI 当外部二进制运行，并通过权限声明限制 `spawn/execute`。
- 平台兼容性符合项目需求：Windows 使用 WebView2，Linux 依赖 WebKitGTK（v2 生态通常要求较新的运行环境）。

## 3. 与 cotton2k 当前代码的契合点

当前仓库已经具备 GUI 化最关键条件：
- Rust 主程序可由 profile 路径直接启动模拟：`/Users/tcztzy/GitHub/cotton2k/src/main.rs:3`。
- 核心入口 `Profile::run()` 已稳定封装：`/Users/tcztzy/GitHub/cotton2k/src/profile.rs:661`。
- 输出约定清晰：结果写入 `output.csv`：`/Users/tcztzy/GitHub/cotton2k/src/profile.rs:944`。
- Python binding 同样是“路径 -> run”调用范式，说明 API 形态简洁，便于映射到 GUI 命令：`/Users/tcztzy/GitHub/cotton2k/bindings/python/src/lib.rs:7`。

结论：**无需改模型算法，即可先做一层 Tauri 应用壳**。

## 4. 方案对比与推荐

### 方案 A：本地 worker（C/S 分进程，每作业一进程，推荐）
- 做法：GUI 作为调度端，调用 `cotton2k-worker run ...`，每个作业启动独立本地进程。
- 优点：
  - 与 GUI 隔离，天然支持 batch 并行，避免进程内全局状态冲突。
  - 失败/取消可通过退出码与 JSONL 事件统一处理。
- 风险：
  - 需要维护 worker 参数协议与打包路径。

### 方案 B：Tauri 进程内直接调用 crate
- 做法：在 `src-tauri` Rust 侧直接调用 `run_job(...)` 或 `Profile::run()`。
- 优点：
  - 调用链最短。
- 风险：
  - 进程内并发运行需要先完成全局状态隔离重构。

### 推荐落地
- **MVP 用方案 A**：worker 分进程 + JSONL 事件流。
- 保留方案 B 作为单作业 fallback。

## 5. 建议的技术架构

### 目录建议
- `apps/desktop/`：前端（Vite + TypeScript，框架可选 React/Vue，建议先 React）。
- `apps/desktop/src-tauri/`：Tauri Rust 后端。
- `src/bin/cotton2k-worker.rs`：本地作业 worker 入口（二进制）。
- `src/bin/cotton2k-batch.rs`：本地批量调度入口（并行拉起多个 worker，可被 GUI 复用）。

### Rust 命令层（Tauri backend）
建议暴露以下命令：
- `validate_profile(path) -> ValidationResult`
- `submit_batch_jobs(requests) -> Vec<JobId>`
- `cancel_job(job_id)`
- `load_output_csv(path) -> TimeseriesPayload`
- `open_output_dir(path)`

并补充事件流：
- `simulation://started`
- `simulation://progress`（按天或阶段）
- `simulation://finished`
- `simulation://failed`

### 前端层（Web UI）
- 参数页：文件选择、基础参数摘要、校验提示。
- 运行页：状态条、日志、取消按钮（先软取消，后续再做强中断）。
- 结果页：多指标折线图 + 表格 + 导出。

### 安全策略（Tauri v2）
- 默认最小权限：仅启用必要文件读取范围。
- 只在确实需要时启用 shell/fs 插件权限，并限制目录 scope。
- 若使用 worker，严格配置 `shell:allow-spawn`/`shell:allow-execute` 的命令和参数白名单。

## 6. 里程碑与工期估算

### M1（3-4 天）：基础骨架
- 初始化 Tauri v2 + 前端模板。
- 接通 `run_simulation` 命令，完成一次端到端运行。
- 验收：可从 GUI 选择 profile 并生成 `output.csv`。

### M2（4-6 天）：可视化与可用性
- 实现 CSV 读取、关键指标图表、错误提示。
- 增加运行元数据（输入文件哈希、开始/结束时间、版本号）。
- 验收：用户可完成“运行 + 看图 + 导出”。

### M3（3-5 天）：打包发布与质量
- 三平台打包脚本（优先 macOS + Windows）。
- 基础回归测试与样例数据 smoke test。
- 验收：生成安装包并能在干净环境启动运行。

> 预估总工期：**2-3 周（1 人）**。

## 7. 风险与应对

- UI 卡顿：
  - 应对：长计算转为 worker 进程执行，前端只做事件订阅与状态管理。
- Linux 构建依赖差异：
  - 应对：在 CI 固定 WebKitGTK 依赖版本并给出安装脚本。
- 结果口径变化引发图表字段不一致：
  - 应对：建立“输出 schema 检查”与字段映射层。
- 自动更新与签名复杂：
  - 应对：MVP 不开自动更新，M3 再引入 updater。

## 8. 验收标准（Definition of Done）

- 用户可在 GUI 中完成一次完整模拟并看到关键曲线。
- 出错时有明确错误提示（含文件路径/字段名）。
- 产物可打包为桌面应用并在目标系统运行。
- 与当前 CLI 在同一输入下，关键输出字段一致（抽样比对）。

## 9. 建议下一步

1. 先在仓库开 `codex/tauri-v2-mvp` 分支。
2. 落地 M1：搭骨架 + 接 `run_simulation`。
3. 我可以继续给你补一份“可直接执行的任务分解清单（Issue 模板 + 验收脚本）”。

---

## 参考资料（调研来源）
- Tauri 2.0 Stable Release (2024-10-02): https://v2.tauri.app/blog/tauri-20/
- Calling Rust from the Frontend（`#[tauri::command]` / async / channels）: https://v2.tauri.app/develop/calling-rust/
- Permissions（v2 权限模型）: https://v2.tauri.app/security/permissions/
- Embedding External Binaries（sidecar / externalBin / shell permission）: https://v2.tauri.app/develop/sidecar/
- Create a Project（v2 初始化方式）: https://v2.tauri.app/start/create-project/
- Webview Versions（平台 WebView 说明）: https://v2.tauri.app/reference/webview-versions/
- Tauri GitHub README（平台与 Linux WebKitGTK 说明）: https://github.com/tauri-apps/tauri
