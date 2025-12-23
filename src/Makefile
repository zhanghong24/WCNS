# 编译器设置
FC = mpif90

# 编译选项说明: 
# -O2: 优化等级2
# -cpp: 启用预处理 (处理 #ifdef 等宏)
# -J build: 将 .mod 模块文件放入 build 目录
# -DPARALLEL: 定义并行宏
# -ffree-line-length-none: 允许代码行无限长 (解决 Line truncated 错误)
# -fdefault-real-8: 将默认实型提升为双精度 (解决 Type mismatch 错误)
# -fdefault-double-8: 将双精度提升为8字节 (与上面配合使用)
FFLAGS = -O2 -cpp -J build -DPARALLEL \
         -ffree-line-length-none \
         -fdefault-real-8 -fdefault-double-8

# 目标文件名称
TARGET = solver
BUILD_DIR = build

# 获取当前目录下的所有 .f90 和 .f 文件
SRCS_F90 = $(wildcard *.f90)
SRCS_F   = $(wildcard *.f)

# 定义目标对象文件列表 (路径前缀加上 build/)
OBJS = $(patsubst %.f90, $(BUILD_DIR)/%.o, $(SRCS_F90)) \
       $(patsubst %.f,   $(BUILD_DIR)/%.o, $(SRCS_F))

# 默认目标
all: $(BUILD_DIR) $(TARGET)

# 创建编译目录
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# 链接规则
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)
	@echo "========================================"
	@echo "编译成功！可执行文件生成于: $(TARGET)"
	@echo "========================================"

# .f90 编译规则
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# .f 编译规则 (针对 mpimod.f 等固定格式文件)
$(BUILD_DIR)/%.o: %.f | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# -----------------------------------------------------------
# 关键模块依赖关系 (确保编译顺序正确)
# -----------------------------------------------------------

# 1. 最底层模块: mpi 包装 和 全局变量
$(BUILD_DIR)/mpimod.o: mpimod.f
$(BUILD_DIR)/Globalvar.o: Globalvar.f90

# 2. 并行模块 (依赖 mpi 和 Globalvar)
$(BUILD_DIR)/parallel.o: parallel.f90 $(BUILD_DIR)/mpimod.o $(BUILD_DIR)/Globalvar.o

# 3. 确保所有其他文件在 Globalvar 和 parallel 编译之后再编译
#    (排除掉基础模块自身，防止循环依赖)
OTHER_OBJS = $(filter-out $(BUILD_DIR)/mpimod.o $(BUILD_DIR)/Globalvar.o $(BUILD_DIR)/parallel.o, $(OBJS))

$(OTHER_OBJS): $(BUILD_DIR)/Globalvar.o $(BUILD_DIR)/parallel.o

# 4. 主程序依赖
$(BUILD_DIR)/Main.o: Main.f90 $(BUILD_DIR)/PreProcess.o $(BUILD_DIR)/RHS.o $(BUILD_DIR)/LHS.o

# 清理规则
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean