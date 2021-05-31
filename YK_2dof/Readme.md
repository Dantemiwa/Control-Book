请参考论文Youla_2_dof设计
以及知乎帖子：https://zhuanlan.zhihu.com/p/376825435/

运行main函数，能看到带滤波器的二自由度控制器的效果。
若不想使用滤波器请将cf_flag变为0
注意由于Wiener-Hopf方法本身的缺陷，sigma不可以设成0，
想取消摄动，请将per_flag设成0
使用了minreal的地方请慎重改动
单自由度的效果可在simulink中查看。
不排除写的函数有问题，欢迎改进指正。