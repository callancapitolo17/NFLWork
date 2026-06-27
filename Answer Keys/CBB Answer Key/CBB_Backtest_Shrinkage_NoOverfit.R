# Targeted anti-overfit evidence: pick n0 from TRAIN ROI, confirm on TEST.
# Plus extreme-subset and EV-bucket breakdowns at the chosen moderate n0.
suppressMessages({ library(duckdb); library(dplyr); library(DBI) })
setwd("/Users/callancapitolo/NFLWork/.claude/worktrees/extreme-samples-shrinkage/Answer Keys")
con <- dbConnect(duckdb(), dbdir="cbb_backtest_shrinkage.duckdb", read_only=TRUE)
res <- dbGetQuery(con, "SELECT * FROM preds"); dbDisconnect(con, shutdown=TRUE)
res$game_date <- as.Date(res$game_date); res <- res[!is.na(res$game_date),]

EV_THRESHOLD<-0.05; BANKROLL<-1000; KELLY_MULT<-0.25
N0_GRID <- c(0,10,20,30,40,50,75,100,150,200,300,500)
calc_ev <- function(p,bo){d<-ifelse(bo>0,1+bo/100,1+100/abs(bo)); p*d-1}
shrink_prob<-function(x,n,pm,n0) if(is.infinite(n0)) pm else (x+n0*pm)/(n+n0)
shrink_var <-function(p,n,n0) if(is.infinite(n0)) rep(0,length(p)) else p*(1-p)/(n+n0+1)
bm_alpha   <-function(e,v){a<-e^2/(e^2+v); ifelse(is.finite(a),a,1)}
kstake     <-function(p,bo,al){d<-ifelse(bo>0,1+bo/100,1+100/abs(bo));b<-d-1;BANKROLL*KELLY_MULT*al*pmax(0,(b*p-(1-p))/b)}

cut_date <- as.Date(quantile(unclass(res$game_date),0.5,type=1), origin="1970-01-01")
train <- res %>% filter(game_date<=cut_date); test <- res %>% filter(game_date>cut_date)

roi_at <- function(df, n0, use_bm=TRUE) {
  p <- shrink_prob(df$x, df$n, df$book_prob, n0)
  ev <- calc_ev(p, df$book_odds); s <- ev>EV_THRESHOLD
  if (sum(s)==0) return(data.frame(n0=n0,n_bets=0,roi=NA,pnl=0))
  d<-df[s,]; pp<-p[s]; dec<-ifelse(d$book_odds>0,1+d$book_odds/100,1+100/abs(d$book_odds))
  al<- if(use_bm) bm_alpha(pp*dec-1, shrink_var(pp,d$n,n0)) else 1
  st<-kstake(pp,d$book_odds,al); pnl<-ifelse(d$actual==1,st*(dec-1),-st)
  data.frame(n0=n0,n_bets=sum(s),roi=round(sum(pnl)/sum(st)*100,2),pnl=round(sum(pnl)))
}

cat("=== n0 chosen on TRAIN by ROI, then CONFIRMED on TEST (no test-peeking) ===\n")
tr <- bind_rows(lapply(N0_GRID, function(n0) roi_at(train,n0)))
te <- bind_rows(lapply(N0_GRID, function(n0) roi_at(test,n0)))
cmp <- data.frame(n0=N0_GRID, train_bets=tr$n_bets, train_roi=tr$roi, train_pnl=tr$pnl,
                  test_bets=te$n_bets, test_roi=te$roi, test_pnl=te$pnl)
print(cmp)
best_train_n0 <- N0_GRID[which.max(tr$pnl)]
cat(sprintf("\nTRAIN-PnL-optimal n0 = %d  -> on TEST: ROI=%.2f%%, PnL=%d (raw TEST: ROI=%.2f%%, PnL=%d)\n\n",
    best_train_n0, te$roi[N0_GRID==best_train_n0], te$pnl[N0_GRID==best_train_n0], te$roi[1], te$pnl[1]))

# Targeted breakdowns at a conservative moderate n0
for (N0 in c(30, 50)) {
  cat(sprintf("================ n0 = %d : targeted breakdowns (TEST) ================\n", N0))
  test2 <- test %>% mutate(p_raw=algo_prob, p_sh=shrink_prob(x,n,book_prob,N0))
  # EV-bucket ROI: raw vs shrunk+BM (the phantom-edge tail lives in high buckets)
  bucket <- function(df,pcol,use_bm){
    d<-df %>% mutate(ev=calc_ev(.data[[pcol]],book_odds)) %>% filter(ev>EV_THRESHOLD)
    if(nrow(d)==0) return(NULL)
    dec<-ifelse(d$book_odds>0,1+d$book_odds/100,1+100/abs(d$book_odds))
    al<- if(use_bm) bm_alpha(d[[pcol]]*dec-1, shrink_var(d[[pcol]],d$n,N0)) else 1
    d$stake<-kstake(d[[pcol]],d$book_odds,al); d$pnl<-ifelse(d$actual==1,d$stake*(dec-1),-d$stake)
    d %>% mutate(evb=cut(ev*100,c(5,10,20,Inf),labels=c("5-10%","10-20%","20%+"))) %>%
      group_by(evb) %>% summarize(n=n(),roi=round(sum(pnl)/sum(stake)*100,1),pnl=round(sum(pnl)),.groups="drop")
  }
  cat("RAW by EV bucket:\n");        print(bucket(test2,"p_raw",FALSE))
  cat("SHRUNK+BM by EV bucket:\n");  print(bucket(test2,"p_sh",TRUE))
  # Extreme subset
  ext <- test2 %>% filter(is_extreme)
  re<-roi_at(ext,0,FALSE); se<-roi_at(ext,N0,TRUE)
  cat(sprintf("EXTREME subset: raw n=%d roi=%s pnl=%s | shrunk+BM n=%d roi=%s pnl=%s\n\n",
      re$n_bets,re$roi,re$pnl, se$n_bets,se$roi,se$pnl))
}
cat("DONE\n")
