# =======================================================================
#  clip_and_align_close.R   ——  对齐 2015-2024 的 Close 价格
# =======================================================================

library(data.table)   # fread 读 csv
library(xts)          # 时间序列对齐
library(lubridate)    # 解析日期

# --------- 1. 配置 -----------------------------------------------------
data_dir   <- "data/stocks"                 # 存 放 CSV 的目录
start_date <- as.Date("2015-01-01")
end_date   <- as.Date("2024-12-31")
out_price_csv <- "output/close_2015_2024.csv"
out_logret_csv <- "output/logret_2015_2024.csv"

# --------- 2. 读取单股票函数 -------------------------------------------
read_one <- function(path) {
  dt <- fread(path, encoding = "UTF-8", showProgress = FALSE)

  # 找到日期列（忽略大小写）
  date_col  <- grep("^date$",   names(dt), ignore.case = TRUE, value = TRUE)[1]
  close_col <- grep("^close$",  names(dt), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(date_col) || is.na(close_col))
    stop("缺少 date 或 close 列: ", basename(path))

  # 解析日期（支持 2002/5/23 或 2015-06-01 等）
  dt[[date_col]] <- gsub("[./]", "-", dt[[date_col]])     # 统一分隔符
  dt[, Date := ymd(dt[[date_col]])]

  # 过滤区间 + 去 NA
  dt <- dt[Date >= start_date & Date <= end_date & !is.na(dt[[close_col]]),
           .(Date, Close = as.numeric(get(close_col)))]

  # 转为 xts，列名设为股票代码（文件名去掉扩展名）
  code <- tools::file_path_sans_ext(basename(path))
  xts(dt$Close, order.by = dt$Date, tzone = "UTC", dimnames = list(NULL, code))
}

# --------- 3. 读取并对齐全部股票 --------------------------------------
csv_files <- list.files(data_dir, "\\.csv$", full.names = TRUE)
stopifnot(length(csv_files) > 0)

price_list <- lapply(csv_files, read_one)

# 取共同交易日（交集）
common_dates <- Reduce(intersect, lapply(price_list, index))
if (length(common_dates) == 0) stop("指定区间没有共同交易日")

price_list <- lapply(price_list, function(x) x[common_dates])
close_xts  <- do.call(merge, c(price_list, join = "inner"))   # m × K

# --------- 4. 计算 Δlog 收益 -------------------------------------------
logret_xts <- diff(log(close_xts))            # 行 = m-1
logret_xts <- na.omit(logret_xts)

# --------- 5. 保存 ------------------------------------------------------
dir.create("output", showWarnings = FALSE)
write.zoo(close_xts, out_price_csv,  sep = ",")
write.zoo(logret_xts, out_logret_csv, sep = ",")

cat("完成：\n",
    "  Close 矩阵 :", nrow(close_xts), "×", ncol(close_xts), " =>", out_price_csv, "\n",
    "  LogRet矩阵 :", nrow(logret_xts), "×", ncol(logret_xts), " =>", out_logret_csv, "\n")

# --------- 6. 若需矩阵喂给 SDE 算法 -----------------------------------
Z_Delta <- coredata(logret_xts)               # 数值矩阵 m × K
ts       <- 1:nrow(Z_Delta) / nrow(Z_Delta)   # 归一化时间

# 示例放进 sim_data 结构后调用
# sim_data <- list(Z_Delta = Z_Delta,
#                  ts      = ts,
#                  n_vec   = rep(1, ncol(Z_Delta)),
#                  mu_true = NULL,
#                  sigma2_true = NULL)
# mu_res <- estimate_mu_from_data(sim_data)
