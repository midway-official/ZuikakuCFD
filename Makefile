# ==========================================================
# 编译器与选项
# ==========================================================

CXX        := mpic++
CXX_STD    := -std=c++17
OPT_FLAGS  := -O3 -march=native -mtune=native -funroll-loops -ffast-math -fomit-frame-pointer
WARN_FLAGS := -Wall -Wextra -Wpedantic -Wshadow -Wno-unused-parameter
CXXFLAGS   := $(CXX_STD) $(OPT_FLAGS) $(WARN_FLAGS)

# ==========================================================
# 目录
# ==========================================================

BUILD_DIR  := build
REPORT_DIR := report

WENO_SRC_DIR := src_WENO
DG_SRC_DIR   := src_DG

# ==========================================================
# 源文件 & 目标文件
# ==========================================================

WENO_SRCS := $(WENO_SRC_DIR)/fluid.cpp $(WENO_SRC_DIR)/solver.cpp
DG_SRCS   := $(DG_SRC_DIR)/fluid.cpp   $(DG_SRC_DIR)/solver.cpp

# 编译产物放在 build/weno/ 和 build/dg/ 下，避免 .o 文件名冲突
WENO_OBJS := $(patsubst $(WENO_SRC_DIR)/%.cpp, $(BUILD_DIR)/weno/%.o, $(WENO_SRCS))
DG_OBJS   := $(patsubst $(DG_SRC_DIR)/%.cpp,   $(BUILD_DIR)/dg/%.o,   $(DG_SRCS))

WENO_DEPS := $(WENO_OBJS:.o=.d)
DG_DEPS   := $(DG_OBJS:.o=.d)

# ==========================================================
# 最终可执行文件
# ==========================================================

TARGET_WENO := solver_WENO
TARGET_DG   := solver_DG

# ==========================================================
# 默认目标：同时构建两个
# ==========================================================

.DEFAULT_GOAL := all

all: $(TARGET_WENO) $(TARGET_DG)
	@echo "✔ 所有目标构建完成：$(TARGET_WENO)  $(TARGET_DG)"

# ==========================================================
# 链接
# ==========================================================

$(TARGET_WENO): $(WENO_OBJS)
	@echo "  链接 $@"
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "✔ $@ 链接成功"

$(TARGET_DG): $(DG_OBJS)
	@echo "  链接 $@"
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "✔ $@ 链接成功"

# ==========================================================
# 编译规则
#
# WENO 源文件 → build/weno/*.o   (头文件从 src_WENO/ 查找)
# DG   源文件 → build/dg/*.o     (头文件从 src_DG/   查找)
# ==========================================================

$(BUILD_DIR)/weno/%.o: $(WENO_SRC_DIR)/%.cpp | $(BUILD_DIR)/weno
	@echo "  编译 [WENO] $<"
	$(CXX) $(CXXFLAGS) -I$(WENO_SRC_DIR) -MMD -MP -c $< -o $@

$(BUILD_DIR)/dg/%.o: $(DG_SRC_DIR)/%.cpp | $(BUILD_DIR)/dg
	@echo "  编译 [DG]   $<"
	$(CXX) $(CXXFLAGS) -I$(DG_SRC_DIR) -MMD -MP -c $< -o $@

# ==========================================================
# 创建目录
# ==========================================================

$(BUILD_DIR)/weno:
	@mkdir -p $@

$(BUILD_DIR)/dg:
	@mkdir -p $@

$(REPORT_DIR):
	@mkdir -p $@

# ==========================================================
# 单独构建目标（方便只编其中一个）
# ==========================================================

weno: $(TARGET_WENO)
	@echo "✔ WENO 构建完成"

dg: $(TARGET_DG)
	@echo "✔ DG 构建完成"

# ==========================================================
# 清理
# ==========================================================

clean:
	@echo "  清理构建产物"
	@rm -rf $(BUILD_DIR) $(TARGET_WENO) $(TARGET_DG)

clean-report:
	@echo "  清理报告"
	@rm -rf $(REPORT_DIR)

distclean: clean clean-report
	@echo "  完全清理"

# ==========================================================
# 自动依赖
# ==========================================================

-include $(WENO_DEPS)
-include $(DG_DEPS)