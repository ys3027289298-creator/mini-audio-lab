# stream

统一输出长度推导，纯头文件，无外部依赖。

| 公式 | 含义 |
|---|---|
| `stream_total_out(N, r) = floor(N/r)` | 整条流的重采样输出总长 |
| `stream_block_out(b, B, r) = floor((b+1)B/r) - floor(bB/r)` | 第 b 块输出，各块之和等于总长 |
| `stream_total_pairs(M) = floor(M/2)` | 整条流的立体声对总数 |
| `stream_block_pairs(n, carry)` | 单块配对数；奇数帧携带到下一块，不补零、不造帧 |

三个处理环节都从这一个头文件取长度，不再各自四舍五入。
