#!/bin/bash

# 批量提交 LSF 作业脚本
# 每个 .lsf 文件提交 3 次

# 设置 .lsf 文件所在目录
LSF_DIR="/home_edu/edu-cmc-sqi22/edu-cmc-sqi22-18/smit/hybrid"

# 定义所有的 LSF 文件
lsf_files=(
    "hybrid_128_4_1.lsf"
    "hybrid_128_4_2.lsf"
    "hybrid_128_4_4.lsf"
    "hybrid_128_4_8.lsf"
    "hybrid_256_8_1.lsf"
    "hybrid_256_8_2.lsf"
    "hybrid_256_8_4.lsf"
    "hybrid_256_8_8.lsf"
    "hybrid_512_8_1.lsf"
    "hybrid_512_8_2.lsf"
    "hybrid_512_8_4.lsf"
    "hybrid_512_8_8.lsf"
)

echo "开始提交 LSF 作业..."
echo "================================"

# 遍历每个 LSF 文件
for lsf_file in "${lsf_files[@]}"; do
    echo ""
    echo "正在处理: $lsf_file"
    echo "--------------------------------"
    
    # 每个文件提交 3 次
    for i in {1..3}; do
        echo "  第 $i 次提交 $lsf_file"
        bsub < "$LSF_DIR/$lsf_file"
        
        # 检查提交是否成功
        if [ $? -eq 0 ]; then
            echo "  ✓ 提交成功"
        else
            echo "  ✗ 提交失败"
        fi
        
        # 稍微延迟，避免过快提交
        sleep 0.5
    done
done

echo ""
echo "================================"
echo "所有作业提交完成！"
echo ""
echo "使用以下命令查看作业状态："
echo "  bjobs        # 查看所有作业"
echo "  bjobs -u all # 查看所有用户的作业"
echo "  bjobs -l     # 查看详细信息"
