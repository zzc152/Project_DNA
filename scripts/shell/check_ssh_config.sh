#!/bin/bash
# 检查 sshd 配置
echo "===== sshd_config 关键配置 ====="
grep -nE "UseDNS|GSSAPIAuth|MaxStartups|LoginGraceTime" /etc/ssh/sshd_config 2>/dev/null
grep -rnE "UseDNS|GSSAPIAuth|MaxStartups|LoginGraceTime" /etc/ssh/sshd_config.d/ 2>/dev/null
echo ""
echo "===== sudo 权限 ====="
sudo -n true 2>&1 && echo "HAS_SUDO" || echo "NO_SUDO"
echo ""
echo "===== 系统信息 ====="
whoami
hostname
uptime
