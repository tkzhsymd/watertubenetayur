real   ,parameter :: alfa=0.05,emin=0.0001  ! iteration check, alfa value is important?
integer,parameter :: itmax=100000        ! iteration max
integer :: ntn(100,10),nb(100),ne(100),ind(100),ntl(100,10) ! data max nodes 100, max tube 100
real    :: ram(100),H(100)                                  ! 

! tube number and node number should be sequential. 1,2,3,4,5....     1,2,4,6... is NG.

! ntn(n,m) 節点nにntn(n,1)個の別の節点が接続されており、その番号がntn(n,2)～ntn(n,ntn(n,1)+1)に記載されている。
! ntl(n,k) 節点nにntl(n,1)本の管路　　が接続されており、その番号がntl(n,2)～ntl(n,ntl(n,1)+1)に記載されている。
! なお、節点ntn(n,2)と節点nを結ぶ管の番号はntl(n,2)となっている。ntn(n,1)=ntl(n,1)に違いない。

! nb(j) j番目の管に接続されている          節点番号
! ne(j) j番目の管に接続されているもう片方の節点番号

! ind(n) n番目の節点が1：計算節点、0：貯水地(const head)なので計算しなくていい

! ram(j) j番目の管の係数　ヘッド差=ram(j)*Q**2 Δh=ram*Q**2

!   H(n) n番目の節点および貯水池のhead, Reservoir(lake) is constant(Boundary condition)


ntn=0; nb=0; ne=0; ram=0.0; nmax=0

 open(10,file='tubes.txt')  ! tube data 管情報（節点情報ではない）
1 continue
 read(10,*,end=10) n,nb(n),ne(n),ram(n) ! node number of both ends of the tube(n), tube constant ram(n)
 nmax=max(nmax,nb(n))
 nmax=max(nmax,ne(n))
 goto 1

10 close(10)
   nl=n     ! 管の本数   total tube number
  non=nmax  ! 節点の総数 total node number

  write(*,*) 'total node number : ',non,' total tube number : ',nl


 do m=1,non ! 節点で掃過する。節点総数は既知か？

 is=1
 do n=1,nl ! 管番号での掃過

 if(nb(n).eq.m)then
 is=is+1
 ntn(m,is)=ne(n)   ! 節点mとつながっている節点番号
 ntl(m,is)=n       ! その管番号
! write(*,*) m,ne(n) ; pause
 endif

 if(ne(n).eq.m)then
 is=is+1
 ntn(m,is)=nb(n)   ! 節点mとつながっている節点番号
 ntl(m,is)=n       ! その管番号
! write(*,*) m,nb(n) ; pause
 endif

 enddo ! 管での掃過

 ntn(m,1)=is-1 ! 節点mはntn(m,1)個の節点とつながっている。
 ntl(m,1)=is-1 ! 何個の管なのか。ntn(m,1)と同じ

! write(*,*) m,ntn(m,1) ; pause

 enddo ! 節点での掃過 all node

 open(10,file='tnumber.dat')
 do n=1,non
!write(*, '(10i5)') n,ntn(n,1),(ntn(n,j),j=2,ntn(n,1)+1) ; pause
 write(10,'(10i5)') n,ntn(n,1),(ntn(n,j),j=2,ntn(n,1)+1)
 write(10,'(10i5)') n,ntl(n,1),(ntl(n,j),j=2,ntl(n,1)+1)
 enddo
 close(10)

 ind=1 ! 計算すべき節点  1 for all nodes but lake node should be 0 because lake water level is constant.

 open(10,file='resor.dat') ! reservoir(lake) data
 read(10,*) nor ! total lake number 貯水地総数（節点のうちのどれか）norはnonの内数
 sum=0.0
 do n=1,nor
 read(10,*) k,H(k); ind(k)=0 ! 節点k番の水位H(k)のデータ、ind=0は与えるから計算しない
 sum=sum+H(k)
 enddo
 ave=sum/nor ! 貯水池の平均水位 ! average lake water level -> initial head of all nodes except lake
 close(10)

 do m=1,non
 if(ind(m).eq.1) H(m)=ave ! 平均水位を全節点(貯水池でなく)の割り当てる。初期値 initial head of each node except lake
 enddo



! caluculation start



! ここまで、節点番号mについて ind(m)=1の水位H(m)について発展計算を行う。ind(m)=0では貯水位であり条件として与える
! その節点mにはntn(m,1)個の別の節点が接続されており、その番号はntn(m,2)...ntn(m,ntn(m,1)+1)に格納されている
! その節点mにはntl(m,1)本の管が接続されており、その番号はntn(m,2)..ntn(m,ntn(m,1)+1)に格納されている。ntn(m,1)=ntl(m,1)
! 全節点数はnonである。
! 節点mと接続されている別の節点（番号ntn(m,n)）との間の管番号はntl(m,n)である。


 do it=1,itmax ! 繰り返し計算する。iteration

 echeck=0.0      ! 繰り返し計算の収束チェッカー convergence checker

 do n=1,non ! 全節点を掃過 all nodes

 sum=0.0    ! 和がゼロにしたい ∑Q initial value

 if(ind(n).eq.1) then ! 1の場合だけ調べる

 do m=1,ntn(n,1) ! nodes connected to node(n)

 i=ntn(n,1+m) ! 相手の節点番号 the node number
 j=ntl(n,1+m) ! その管番号     the tube number
 dh=H(i)-H(n) ! この値が負のときが基本 head difference between i,n
 Q=sqrt(abs(dh)/ram(j)) ! absolute value of flow quantitity no direction

! if conceptual Q direction is from n->i. if i->n ( H(i)>H(n);dh>0 ), then -Q

  if(dh.gt.0.0) Q=-Q ! Qは出ていく方(n->i)が正ときめたのでH(n)>H(i)が基本。その逆

 sum=sum+Q    ! 足していく suming 

 enddo ! 節点nに関係している別の節点・管路で掃過

 echeck=amax1(abs(sum),echeck) ! sumはゼロであって欲しいもの その最大値を調べる

! if(sum>0) too much going out, H(n) goes down
! if(sum<0) too much coming in, H(n) goew up

 H(n)=H(n)-alfa*sum ! sumの値に応じて変更する。sum>0なら減らし、sum<0なら増やす。

 endif

 enddo ! 全節点 for all nodes

if(echeck.lt.emin) then ! 収束、計算終了、最終処理　水位と流量記録

! just output the results

open(10,file='output.dat')

do m=1,nl ! 管ごとに出力
 h1=H(nb(m)); h2=H(ne(m)); dh=h1-h2
 Q=sqrt(abs(dh)/ram(m))  
 if(dh.lt.0.0) Q=-Q      
 write(10,'(3i5,2f8.2,f9.4)') m, nb(m), ne(m), h1, h2, Q
enddo

close(10)

stop

endif

enddo ! 収束繰り返し計算 iteration

print *, 'no convergence'

open(10,file='outpute.dat')

do m=1,nl ! 管ごとに出力
 h1=H(nb(m)); h2=H(ne(m)); dh=h1-h2
 Q=sqrt(abs(dh)/ram(m))  
 if(dh.lt.0.0) Q=-Q      
 write(10,'(3i5,2f8.2,f9.4)') m, nb(m), ne(m), h1, h2, Q
enddo

close(10)

stop
end
