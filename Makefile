# ==========================================================
# 编译器与选项
# ==========================================================

CXX        := mpic++
CXX_STD     := -std=c++17
OPT_FLAGS   := -O3 -march=native -mtune=native -funroll-loops 
WARN_FLAGS  := -Wall -Wextra -Wpedantic -Wshadow -Wno-unused-parameter
CXXFLAGS    := $(CXX_STD) $(OPT_FLAGS) $(WARN_FLAGS)

# ==========================================================
# 目录
# ==========================================================

BUILD_DIR   := build
REPORT_DIR  := report

WENO_SRC_DIR := src_WENO
DG_SRC_DIR   := src_DG

# ==========================================================
# 源文件 & 宏参数定义
# ==========================================================

WENO_SRCS   := $(WENO_SRC_DIR)/fluid.cpp $(WENO_SRC_DIR)/solver.cpp
DG_SRCS     := $(DG_SRC_DIR)/fluid.cpp   $(DG_SRC_DIR)/solver.cpp

# 定义传递给代码的宏 DG_P_VAL 的取值范围 (P=1, 2, 3)
DG_P_VALUES := 1 2 3

# ==========================================================
# 最终可执行文件
# ==========================================================

TARGET_WENO := solver_WENO
# 动态生成目标名：solver_DG2 solver_DG3 solver_DG4
# 计算公式：后缀 = DG_P_VAL + 1
TARGETS_DG  := $(foreach p, $(DG_P_VALUES), solver_DG$(shell echo $$(($(p)+1))))

# ==========================================================
# 默认目标
# ==========================================================

.DEFAULT_GOAL := all

all: $(TARGET_WENO) $(TARGETS_DG)
	@echo "✔ 所有目标构建完成：$(TARGET_WENO) $(TARGETS_DG)"

# ==========================================================
# WENO 编译规则
# ==========================================================

WENO_OBJS := $(patsubst $(WENO_SRC_DIR)/%.cpp, $(BUILD_DIR)/weno/%.o, $(WENO_SRCS))
-include $(WENO_OBJS:.o=.d)

$(TARGET_WENO): $(WENO_OBJS)
	@echo "  链接 $@"
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "✔ $@ 链接成功"

$(BUILD_DIR)/weno/%.o: $(WENO_SRC_DIR)/%.cpp | $(BUILD_DIR)/weno
	@echo "  编译 [WENO] $<"
	$(CXX) $(CXXFLAGS) -I$(WENO_SRC_DIR) -MMD -MP -c $< -o $@

# ==========================================================
# DG 编译规则 (自动化生成 solver_DG2, DG3, DG4)
# ==========================================================

# 函数模板
# $(1): 传入代码的 DG_P_VAL (1, 2, 或 3)
# $(2): 对应的阶数后缀 (2, 3, 或 4)
define DG_RULE_TEMPLATE
DG$(2)_OBJS := $$(patsubst $$(DG_SRC_DIR)/%.cpp, $$(BUILD_DIR)/dg$(2)/%.o, $$(DG_SRCS))

-include $$(DG$(2)_OBJS:.o=.d)

solver_DG$(2): $$(DG$(2)_OBJS)
	@echo "  链接 $$@"
	$$(CXX) $$(CXXFLAGS) $$^ -o $$@
	@echo "✔ $$@ 链接成功"

$$(BUILD_DIR)/dg$(2)/%.o: $$(DG_SRC_DIR)/%.cpp | $$(BUILD_DIR)/dg$(2)
	@echo "  编译 [DG P=$(1), Order=$(2)] $$<"
	$$(CXX) $$(CXXFLAGS) -DDG_P_VAL=$(1) -I$$(DG_SRC_DIR) -MMD -MP -c $$< -o $$@

$$(BUILD_DIR)/dg$(2):
	@mkdir -p $$@
endef

# 实例化规则：循环传入 P 及其对应的阶数 (P+1)
$(foreach p, $(DG_P_VALUES), $(eval $(call DG_RULE_TEMPLATE,$(p),$(shell echo $$(($(p)+1))))))

# ==========================================================
# 通用目录与清理
# ==========================================================

$(BUILD_DIR)/weno:
	@mkdir -p $@

$(REPORT_DIR):
	@mkdir -p $@

weno: $(TARGET_WENO)
dg: $(TARGETS_DG)

clean:
	@echo "  清理构建产物"
	@rm -rf $(BUILD_DIR) $(TARGET_WENO) $(TARGETS_DG)

clean-report:
	@echo "  清理报告"
	@rm -rf $(REPORT_DIR)

distclean: clean clean-report
	@echo "  完全清理"

.PHONY: all weno dg clean clean-report distclean