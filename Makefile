# ==========================================================
# 编译器与选项
# ==========================================================

CXX      := mpic++
CXX_STD  := -std=c++17
OPT_FLAGS := -O3 -march=native -mtune=native -funroll-loops -ffast-math -fomit-frame-pointer
WARN_FLAGS := -Wall -Wextra -Wpedantic -Wshadow -Wno-unused-parameter
CXXFLAGS := $(CXX_STD) $(OPT_FLAGS) $(WARN_FLAGS)
INCLUDES := -Isrc

# ==========================================================
# 目录
# ==========================================================

SRC_DIR   := src
BUILD_DIR := build
REPORT_DIR := report

# ==========================================================
# 目标程序
# ==========================================================

TARGET    := solver
SRCS      := $(SRC_DIR)/fluid.cpp $(SRC_DIR)/solver.cpp
OBJS      := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS      := $(OBJS:.o=.d)

# ==========================================================
# 默认目标
# ==========================================================

.DEFAULT_GOAL := all

all: $(BUILD_DIR) $(TARGET)
	@echo " 所有目标构建完成"

$(TARGET): $(OBJS)
	@echo " 链接 $@"
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo " $@ 链接成功"

# ==========================================================
# 编译规则（自动依赖）
# ==========================================================

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@echo " 编译 $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

# ==========================================================
# 创建目录
# ==========================================================

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(REPORT_DIR)

# ==========================================================
# 清理
# ==========================================================

clean:
	@echo " 清理构建产物"
	@rm -rf $(BUILD_DIR) $(TARGET)

clean-report:
	@echo " 清理报告"
	@rm -rf $(REPORT_DIR)

distclean: clean clean-report
	@echo " 完全清理"

# ==========================================================
# 自动依赖
# ==========================================================

-include $(DEPS)