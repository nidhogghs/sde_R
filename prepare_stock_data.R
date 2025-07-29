# =========================================================
#  prepare_stock_data.R  —— 批量裁切 & 对齐 & 导出 CSV
# =========================================================

library(data.table)   # 高效读写
library(xts)          # 时间序列对齐
library(lubridate)    # 日期解析

# ---------- 配置区域 ----------------------------------------------------------
data_dir   <- "data/stocks"            # 存放原始 CSV 的文件夹
price_col  <- "close"                  # 收盘价列名（忽略大小写）
start_date <- as.Date("2015-01-01")    # 起始
end_date   <- as.Date("2024-12-31")    # 结束
min_obs    <- 200                      # 至少有多少交易日才保留
out_dir    <- "output"                 # 对齐后 CSV 输出目录
close_csv  <- file.path(out_dir, "close_2015_2024.csv")
logret_csv <- file.path(out_dir, "logret_2015_2024.csv")
# -----------------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

csv_files <- list.files(data_dir, "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("未在 ", data_dir, " 找到任何 CSV 文件")

# ---- 单股票读取函数 ----------------------------------------------------------
read_one <- function(path) {
  fn <- basename(path)
  dt <- tryCatch(
    fread(path, select = c("date", price_col), showProgress = FALSE),
    error = function(e) {
      message("⚠ 无法读取，跳过：", fn)
      return(NULL)
    }
  )
  if (is.null(dt)) return(NULL)

  setnames(dt, tolower(names(dt)))              # -> date / close
  dt[, date := ymd(gsub("[./]", "-", date))]    # 兼容 2002/5/23、2002-05-23
  dt <- dt[date >= start_date & date <= end_date & !is.na(close),
           .(date, close = as.numeric(close))]

  if (nrow(dt) < min_obs) {
    message("⚠ 记录不足 ", min_obs, " 行，跳过：", fn)
    return(NULL)
  }

  code <- tools::file_path_sans_ext(fn)
  xts(dt$close, order.by = dt$date, tzone = "UTC", dimnames = list(NULL, code))
}

# ---- 批量读取 ---------------------------------------------------------------
message("📥 读取 ", length(csv_files), " 只股票…")
price_list <- lapply(csv_files, read_one)
price_list <- Filter(Negate(is.null), price_list)
if (length(price_list) < 2) stop("有效股票少于 2 只，请检查数据")
message("✔ 有效股票数：", length(price_list))

# ---- 交集对齐 ---------------------------------------------------------------
merge_in <- function(x, y) merge(x, y, join = "inner")
close_xts <- Reduce(merge_in, price_list)   # m × K
if (NCOL(close_xts) < 2) stop("交集后仍不足 2 只股票")
message("✔ 共同交易日数：", nrow(close_xts))

# ---- Δlog 收益 --------------------------------------------------------------
logret_xts <- diff(log(close_xts))[-1]      # 去掉首行 NA

# ---- 保存 --------------------------------------------------------------------
write.zoo(close_xts, close_csv, sep = ",")
write.zoo(logret_xts, logret_csv, sep = ",")

cat("已写出：\n",
    "   • 收盘价矩阵：", close_csv, "\n",
    "   • Δlog 收益  ：", logret_csv, "\n",
    "   维度：Close=", nrow(close_xts), "×", ncol(close_xts),
    "；LogRet=", nrow(logret_xts), "×", ncol(logret_xts), "\n")
