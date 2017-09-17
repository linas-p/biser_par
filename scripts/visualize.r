p_list <- c(1, 2, 3, 4, 5, 6, 7, 8);
N_list <- c(16, 32, 48, 64, 96, 128, 192, 256);

######
dat <- c()

for(p in p_list){
  for(n in N_list){
    dat <- rbind(dat, c(p, n, read.table(paste("../output/N_", n, "_p_", p, sep = ""))$V1));
  }
}
dat <- as.data.frame(dat);
names(dat) <- c("n", "N", "Laikas");




######
dat_mpi <- c()

for(p in p_list){
  for(n in N_list){
    dat_mpi <- rbind(dat_mpi, c(p, n, read.table(paste("../output/mpi_N_", n, "_p_", p, sep = ""))$V1));
  }
}
dat_mpi <- as.data.frame(dat_mpi);
names(dat_mpi) <- c("n", "N", "Laikas");


dat_omp <- c()

for(p in p_list){
  for(n in N_list){
    dat_omp <- rbind(dat_omp, c(p, n, read.table(paste("../output/omp_N_", n, "_p_", p, sep = ""))$V1));
  }
}
dat_omp <- as.data.frame(dat_omp);
names(dat_omp) <- c("n", "N", "Laikas");


dat_ref <- c()

for(n in N_list){
  dat_ref <- rbind(dat_ref, c(1, n, read.table(paste("../output/ref_N_", n, "_p_", 1, sep = ""))$V1));
}

dat_ref <- as.data.frame(dat_ref);
names(dat_ref) <- c("n", "N", "Laikas");







n_n <- length(N_list)
nn <- 1:n_n
plot(dat_mpi[nn,2], dat_mpi[nn,3], log = 'x', type = 'l', ylim = c(0, 10), xlab = "Mazgu dydis n", ylab = "laikas (ms)", 
     main = "Srendimo laiko priklausomybė")
lines(dat_mpi[nn+n_n,2], dat_mpi[nn+n_n,3], col = 'red')
lines(dat_mpi[nn+n_n*2,2], dat_mpi[nn+n_n*2,3], col = 'blue')
lines(dat_mpi[nn+n_n*3,2], dat_mpi[nn+n_n*3,3], col = 'green')


setEPS()
postscript("steady.eps")
co <- topo.colors(8);
plot(dat$n[dat$N == 16], (p_list*dat$Laikas[dat$N == 16][1])/dat$Laikas[dat$N == 16], log = 'xy', type = 'o', 
     ylim = c(1, 8), xlab = "procesoriu skaicius #n", ylab = TeX("Spartinimo koeficientas $S = T_1/T_N$"), col = co[1], 
     main = "Spartinimo koeficientas")
lines(dat$n[dat$N == 32], (p_list*dat$Laikas[dat$N == 32][1])/dat$Laikas[dat$N == 32], col = co[2])
lines(dat$n[dat$N == 48], (p_list*dat$Laikas[dat$N == 48][1])/dat$Laikas[dat$N == 48], col = co[3])
lines(dat$n[dat$N == 64], (p_list*dat$Laikas[dat$N == 64][1])/dat$Laikas[dat$N == 64], col = co[4])
lines(dat$n[dat$N == 96], (p_list*dat$Laikas[dat$N == 96][1])/dat$Laikas[dat$N == 96], col = co[5])
lines(dat$n[dat$N == 128 ], (p_list*dat$Laikas[dat$N == 128 ][1])/dat$Laikas[dat$N == 128 ], col = co[6])
lines(dat$n[dat$N == 192], (p_list*dat$Laikas[dat$N == 192][1])/dat$Laikas[dat$N == 192], col = co[7])
lines(dat$n[dat$N == 256], (p_list*dat$Laikas[dat$N == 256][1])/dat$Laikas[dat$N == 256], col = co[8])
lines(p_list, p_list, type = 'o')
legend("topleft", c("N = 16", "N = 32", "N = 48", "N = 64", "N = 96", "N = 128", "N = 192", "N = 256", "Idialus"), 
       lty = rep(1, 9), col = c(co, "black"));
dev.off()


setEPS()
postscript("e.eps")
co <- topo.colors(8);
plot(dat$n[dat$N == 16], (dat$Laikas[dat$N == 16][1])/dat$Laikas[dat$N == 16], log = 'x', type = 'o', 
     ylim = c(0, 1.0), xlab = "procesoriu skaicius #n", ylab = TeX("Efektyvumas"), col = co[1], 
     main = "Efektyvumas (MPI)")
lines(dat$n[dat$N == 32], (dat$Laikas[dat$N == 32][1])/dat$Laikas[dat$N == 32], col = co[2])
lines(dat$n[dat$N == 48], (dat$Laikas[dat$N == 48][1])/dat$Laikas[dat$N == 48], col = co[3])
lines(dat$n[dat$N == 64], (dat$Laikas[dat$N == 64][1])/dat$Laikas[dat$N == 64], col = co[4])
lines(dat$n[dat$N == 96], (dat$Laikas[dat$N == 96][1])/dat$Laikas[dat$N == 96], col = co[5])
lines(dat$n[dat$N == 128 ], (dat$Laikas[dat$N == 128 ][1])/dat$Laikas[dat$N == 128 ], col = co[6])
lines(dat$n[dat$N == 192], (dat$Laikas[dat$N == 192][1])/dat$Laikas[dat$N == 192], col = co[7])
lines(dat$n[dat$N == 256], (dat$Laikas[dat$N == 256][1])/dat$Laikas[dat$N == 256], col = co[8])
lines(c(1, 2, 4, 8), c(1, 2, 4, 8)/(c(1, 2, 4, 8)))
legend("bottomleft", c("N = 16", "N = 32", "N = 48", "N = 64", "N = 96", "N = 128", "N = 192", "N = 256", "Idialus"), 
       lty = rep(1, 9), col = c(co, "black"));

dev.off()







setEPS()
postscript("s_mpi.eps")
co <- topo.colors(8);
plot(dat_mpi$n[dat_mpi$N == 16], (p_list*dat_ref$Laikas[dat_ref$N == 16][1])/dat_mpi$Laikas[dat_mpi$N == 16], log = 'xy', type = 'o', 
     ylim = c(1, 8), xlab = "procesoriu skaicius #n", ylab = TeX("Spartinimo koeficientas $S = T_1/T_N$"), col = co[1], 
     main = "Spartinimo koeficientas(MPI)")
lines(dat_mpi$n[dat_mpi$N == 32], (p_list*dat_ref$Laikas[dat_ref$N == 32][1])/dat_mpi$Laikas[dat_mpi$N == 32], col = co[2])
lines(dat_mpi$n[dat_mpi$N == 48], (p_list*dat_ref$Laikas[dat_ref$N == 48][1])/dat_mpi$Laikas[dat_mpi$N == 48], col = co[3])
lines(dat_mpi$n[dat_mpi$N == 64], (p_list*dat_ref$Laikas[dat_ref$N == 64][1])/dat_mpi$Laikas[dat_mpi$N == 64], col = co[4])
lines(dat_mpi$n[dat_mpi$N == 96], (p_list*dat_ref$Laikas[dat_ref$N == 96][1])/dat_mpi$Laikas[dat_mpi$N == 96], col = co[5])
lines(dat_mpi$n[dat_mpi$N == 128 ], (p_list*dat_ref$Laikas[dat_ref$N == 128 ][1])/dat_mpi$Laikas[dat_mpi$N == 128 ], col = co[6])
lines(dat_mpi$n[dat_mpi$N == 192], (p_list*dat_ref$Laikas[dat_ref$N == 192][1])/dat_mpi$Laikas[dat_mpi$N == 192], col = co[7])
lines(dat_mpi$n[dat_mpi$N == 256], (p_list*dat_ref$Laikas[dat_ref$N == 256][1])/dat_mpi$Laikas[dat_mpi$N == 256], col = co[8])
lines(p_list, p_list, type = 'o')
legend("topleft", c("N = 16", "N = 32", "N = 48", "N = 64", "N = 96", "N = 128", "N = 192", "N = 256", "Idialus"), 
       lty = rep(1, 9), col = c(co, "black"));
dev.off()

setEPS()
postscript("e_mpi.eps")
co <- topo.colors(8);
plot(dat_mpi$n[dat_mpi$N == 16], (dat_ref$Laikas[dat_ref$N == 16][1])/dat_mpi$Laikas[dat_mpi$N == 16], log = 'x', type = 'o', 
     ylim = c(0, 1.0), xlab = "procesoriu skaicius #n", ylab = TeX("Efektyvumas"), col = co[1], 
     main = "Efektyvumas (MPI)")
lines(dat_mpi$n[dat_mpi$N == 32], (dat_ref$Laikas[dat_ref$N == 32][1])/dat_mpi$Laikas[dat_mpi$N == 32], col = co[2])
lines(dat_mpi$n[dat_mpi$N == 48], (dat_ref$Laikas[dat_ref$N == 48][1])/dat_mpi$Laikas[dat_mpi$N == 48], col = co[3])
lines(dat_mpi$n[dat_mpi$N == 64], (dat_ref$Laikas[dat_ref$N == 64][1])/dat_mpi$Laikas[dat_mpi$N == 64], col = co[4])
lines(dat_mpi$n[dat_mpi$N == 96], (dat_ref$Laikas[dat_ref$N == 96][1])/dat_mpi$Laikas[dat_mpi$N == 96], col = co[5])
lines(dat_mpi$n[dat_mpi$N == 128 ], (dat_ref$Laikas[dat_ref$N == 128 ][1])/dat_mpi$Laikas[dat_mpi$N == 128 ], col = co[6])
lines(dat_mpi$n[dat_mpi$N == 192], (dat_ref$Laikas[dat_ref$N == 192][1])/dat_mpi$Laikas[dat_mpi$N == 192], col = co[7])
lines(dat_mpi$n[dat_mpi$N == 256], (dat_ref$Laikas[dat_ref$N == 256][1])/dat_mpi$Laikas[dat_mpi$N == 256], col = co[8])
lines(c(1, 2, 4, 8), c(1, 2, 4, 8)/(c(1, 2, 4, 8)))
legend("bottomleft", c("N = 16", "N = 32", "N = 48", "N = 64", "N = 96", "N = 128", "N = 192", "N = 256", "Idialus"), 
       lty = rep(1, 9), col = c(co, "black"));

dev.off()

setEPS()
postscript("s_omp.eps")
co <- topo.colors(8);
plot(dat_omp$n[dat_omp$N == 16], (p_list*dat_ref$Laikas[dat_ref$N == 16][1])/dat_omp$Laikas[dat_omp$N == 16], log = 'xy', type = 'o', 
     ylim = c(1, 8), xlab = "procesoriu skaicius #n", ylab = TeX("Spartinimo koeficientas $S = T_1/T_N$"), col = co[1], 
     main = "Spartinimo koeficientas (OMP)")
lines(dat_omp$n[dat_omp$N == 32], (p_list*dat_ref$Laikas[dat_ref$N == 32][1])/dat_omp$Laikas[dat_omp$N == 32], col = co[2])
lines(dat_omp$n[dat_omp$N == 48], (p_list*dat_ref$Laikas[dat_ref$N == 48][1])/dat_omp$Laikas[dat_omp$N == 48], col = co[3])
lines(dat_omp$n[dat_omp$N == 64], (p_list*dat_ref$Laikas[dat_ref$N == 64][1])/dat_omp$Laikas[dat_omp$N == 64], col = co[4])
lines(dat_omp$n[dat_omp$N == 96], (p_list*dat_ref$Laikas[dat_ref$N == 96][1])/dat_omp$Laikas[dat_omp$N == 96], col = co[5])
lines(dat_omp$n[dat_omp$N == 128 ], (p_list*dat_ref$Laikas[dat_ref$N == 128 ][1])/dat_omp$Laikas[dat_omp$N == 128 ], col = co[6])
lines(dat_omp$n[dat_omp$N == 192], (p_list*dat_ref$Laikas[dat_ref$N == 192][1])/dat_omp$Laikas[dat_omp$N == 192], col = co[7])
lines(dat_omp$n[dat_omp$N == 256], (p_list*dat_ref$Laikas[dat_ref$N == 256][1])/dat_omp$Laikas[dat_omp$N == 256], col = co[8])
lines(1:8, 1:8, type = 'o');
legend("topleft", c("N = 16", "N = 32", "N = 48", "N = 64", "N = 96", "N = 128", "N = 192", "N = 256", "Idialus"), 
       lty = rep(1, 9), col = c(co, "black"));
dev.off()

setEPS()
postscript("e_omp.eps")
co <- topo.colors(8);
plot(dat_omp$n[dat_omp$N == 16], (dat_ref$Laikas[dat_ref$N == 16][1])/dat_omp$Laikas[dat_omp$N == 16], log = 'x', type = 'o', 
     ylim = c(0, 1.0), xlab = "procesoriu skaicius #n", ylab = TeX("Efektyvumas"), col = co[1], 
     main = "Efektyvumas (OMP)")
lines(dat_omp$n[dat_omp$N == 32], (dat_ref$Laikas[dat_ref$N == 32][1])/dat_omp$Laikas[dat_omp$N == 32], col = co[2])
lines(dat_omp$n[dat_omp$N == 48], (dat_ref$Laikas[dat_ref$N == 48][1])/dat_omp$Laikas[dat_omp$N == 48], col = co[3])
lines(dat_omp$n[dat_omp$N == 64], (dat_ref$Laikas[dat_ref$N == 64][1])/dat_omp$Laikas[dat_omp$N == 64], col = co[4])
lines(dat_omp$n[dat_omp$N == 96], (dat_ref$Laikas[dat_ref$N == 96][1])/dat_omp$Laikas[dat_omp$N == 96], col = co[5])
lines(dat_omp$n[dat_omp$N == 128 ], (dat_ref$Laikas[dat_ref$N == 128 ][1])/dat_omp$Laikas[dat_omp$N == 128 ], col = co[6])
lines(dat_omp$n[dat_omp$N == 192], (dat_ref$Laikas[dat_ref$N == 192][1])/dat_omp$Laikas[dat_omp$N == 192], col = co[7])
lines(dat_omp$n[dat_omp$N == 256], (dat_ref$Laikas[dat_ref$N == 256][1])/dat_omp$Laikas[dat_omp$N == 256], col = co[8])
lines(1:8, 1:8/1:8, type = 'o');
legend("bottomleft", c("N = 16", "N = 32", "N = 48", "N = 64", "N = 96", "N = 128", "N = 192", "N = 256", "Idialus"), 
       lty = rep(1, 9), col = c(co, "black"));
dev.off()

