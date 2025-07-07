# OFDM接收机模块化系统

本文件夹包含了一个完整的模块化OFDM接收机系统，将原来的单体函数重构为多个独立的功能模块。

## 文件结构

### 核心模块文件
- `ofdm_receiver_main.m` - 主控制模块，协调各个子模块的执行
- `ofdm_receiver_init.m` - 接收机初始化模块，处理输入信号预处理
- `ofdm_generate_local_signals.m` - 本地参考信号生成模块
- `ofdm_tracking_init.m` - 跟踪环路初始化模块，设置DLL和PLL参数
- `ofdm_signal_tracking.m` - 信号跟踪处理模块，执行码跟踪和载波跟踪
- `ofdm_results_output.m` - 结果输出和可视化模块

### 测试和示例
- `test_ofdm_receiver.m` - 完整的测试脚本，演示如何使用模块化系统

## 快速开始

### 1. 查看系统架构图

首先建议查看系统架构图来理解整个系统的结构：

```bash
# 方法1：在浏览器中打开（推荐）
open system_block_diagram.svg

# 方法2：在VSCode中查看
# 直接点击VSCode中的 system_block_diagram.svg 标签页

# 方法3：在Finder中查看（macOS）
open .  # 打开当前文件夹，然后双击SVG文件
```

配套的详细说明文档：`SYSTEM_ARCHITECTURE.md`

### 2. 基本使用方法

```matlab
% 运行测试脚本
cd ofdm_receiver_modules
test_ofdm_receiver
```

### 2. 自定义使用

```matlab
% 1. 初始化仿真参数
simSettings = init(-15);  % SNR = -15 dB
simSettings.code = your_code_sequence;

% 2. 接收机初始化
[processed_signal, init_params] = ofdm_receiver_init(simSettings, received_signal);

% 3. 生成本地参考信号
local_signals = ofdm_generate_local_signals(simSettings);

% 4. 初始化跟踪参数
tracking_params = ofdm_tracking_init(simSettings);

% 5. 执行信号跟踪
for loop_idx = 1:num_loops
    tracking_result = ofdm_signal_tracking(processed_signal, local_signals, ...
                                          tracking_params, loop_idx);
    all_results{loop_idx} = tracking_result;
end

% 6. 输出结果
output_options.save_figures = true;
output_options.show_plots = true;
ofdm_results_output(all_results, simSettings, output_options);
```

## 模块详细说明

### ofdm_receiver_main.m
**功能**: 主控制模块
- 协调各个子模块的执行顺序
- 管理数据流和状态传递
- 提供统一的系统接口

**输入参数**:
- `simSettings`: 仿真设置参数
- `received_signal`: 接收到的信号

**输出参数**:
- `final_results`: 完整的跟踪结果

### ofdm_receiver_init.m
**功能**: 接收机初始化
- 输入参数验证和信号预处理
- 信号格式识别和转换（支持I/Q分离和复数格式）
- 信号长度检查和调整

**关键特性**:
- 自动识别信号格式
- 信号长度自适应调整
- 信号质量初步检查

### ofdm_generate_local_signals.m
**功能**: 本地参考信号生成
- 生成OFDM调制的本地参考信号
- 信号插值和保护间隔处理
- 多频带信号支持

**关键特性**:
- 支持多种码类型（双极性、单极性等）
- 自动插值和保护间隔添加
- 信号质量统计和PAPR计算

### ofdm_tracking_init.m
**功能**: 跟踪环路初始化
- DLL（码跟踪环路）参数设置
- PLL（载波跟踪环路）参数设置
- 环路滤波器系数计算

**关键特性**:
- 根据信噪比自适应调整环路带宽
- 支持FLL辅助的PLL结构
- 完整的参数合理性检查

### ofdm_signal_tracking.m
**功能**: 信号跟踪处理
- 码跟踪（早晚门鉴相器）
- 载波跟踪（Costas环路/FLL）
- 实时性能监控

**关键特性**:
- 多频带相关处理
- FLL到PLL的自动切换
- 实时SNR估计和失锁检测

### ofdm_results_output.m
**功能**: 结果输出和可视化
- 跟踪性能分析和统计
- 多种可视化图表生成
- 数据和报告文件保存

**输出内容**:
- 跟踪误差时间序列图
- 频率跟踪图
- 相关器输出分析图
- 性能监控图
- 统计分析图
- 详细的文本报告

## 系统优势

1. **模块化设计**: 每个模块职责明确，易于理解和维护
2. **标准化接口**: 模块间通过标准化的数据结构通信
3. **独立测试**: 每个模块可以独立测试和调试
4. **易于扩展**: 可以方便地添加新功能或替换特定模块
5. **详细监控**: 提供丰富的处理过程信息和性能指标
6. **可视化支持**: 内置完整的结果分析和可视化功能

## 依赖要求

### 必需的函数文件
确保以下函数在MATLAB路径中可用：
- `init.m` - 初始化函数
- `generateCrossOFDM.m` - OFDM调制函数
- `interpo.m` - 插值函数
- `calLoopCoef.m` - 环路系数计算函数
- `weil10230_signed.mat` - 本地码数据文件

### MATLAB工具箱
- Signal Processing Toolbox（用于信号处理函数）
- Statistics and Machine Learning Toolbox（用于统计分析）

## 配置选项

### 跟踪参数调整
可以通过修改`ofdm_tracking_init.m`中的参数来调整跟踪性能：
- DLL噪声带宽：影响码跟踪精度
- PLL噪声带宽：影响载波跟踪精度
- 阻尼比：影响环路稳定性
- FLL切换阈值：控制FLL到PLL的切换时机

### 输出选项配置
在调用`ofdm_results_output.m`时可以设置：
```matlab
output_options.save_figures = true;    % 保存图形文件
output_options.show_plots = true;      % 显示图形
output_options.save_data = true;       % 保存数据文件
output_options.output_dir = './results'; % 输出目录
output_options.file_prefix = 'ofdm';   % 文件名前缀
```

## 故障排除

### 常见问题
1. **"函数未找到"错误**: 确保所有依赖函数在MATLAB路径中
2. **内存不足**: 对于长信号，考虑分段处理
3. **跟踪性能差**: 调整环路带宽参数或检查信号质量
4. **文件路径错误**: 确认 `weil10230_signed.mat` 文件路径正确

### 调试建议
1. 使用`test_fixed_system.m`验证修复后的系统功能
2. 使用`test_ofdm_receiver.m`进行详细的模块测试
3. 检查每个模块的输出信息和警告
4. 使用可视化结果分析跟踪性能

### 系统修复状态 (2025年7月)
✅ **已修复的问题**:
1. 添加了缺失的参数字段 (`fi`, `fp`, `dt`, `t_total`)
2. 创建了 `ofdm_demodulation.m` 模块
3. 创建了 `ofdm_performance_evaluation.m` 模块  
4. 创建了 `ofdm_visualization.m` 模块
5. 修复了主函数中的计时器问题
6. 修复了MATLAB语法兼容性问题
7. **修复了信号生成中的维度不匹配问题**:
   - 修复了 `generateCrossOFDM.m` 中 `time_data` 转置导致的维度错误
   - 添加了 `reshape` 操作的安全检查和错误处理
   - 改进了信号连接时的维度一致性处理
   - 增强了插值函数的维度处理能力

✅ **系统状态**: 所有核心模块已完整实现并通过测试

### 测试脚本
- `test_fixed_system.m` - 完整系统功能测试
- `test_dimension_fix.m` - 专门测试维度修复的脚本
- `test_ofdm_receiver.m` - 原始模块测试脚本

## 版本信息
- 版本: 1.0
- 创建日期: 2025年7月
- 作者: OFDM接收机开发

## 许可证
本代码仅供学习和研究使用。