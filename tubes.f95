real   ,parameter :: alfa=0.05,emin=0.0001  ! iteration check, alfa value is important?
integer,parameter :: itmax=100000        ! iteration max
integer :: ntn(100,10),nb(100),ne(100),ind(100),ntl(100,10) ! data max nodes 100, max tube 100
real    :: ram(100),H(100)                                  ! 

! tube number and node number should be sequential. 1,2,3,4,5....     1,2,4,6... is NG.

! ntn(n,m) �ߓ_n��ntn(n,1)�̕ʂ̐ߓ_���ڑ�����Ă���A���̔ԍ���ntn(n,2)�`ntn(n,ntn(n,1)+1)�ɋL�ڂ���Ă���B
! ntl(n,k) �ߓ_n��ntl(n,1)�{�̊ǘH�@�@���ڑ�����Ă���A���̔ԍ���ntl(n,2)�`ntl(n,ntl(n,1)+1)�ɋL�ڂ���Ă���B
! �Ȃ��A�ߓ_ntn(n,2)�Ɛߓ_n�����Ԋǂ̔ԍ���ntl(n,2)�ƂȂ��Ă���Bntn(n,1)=ntl(n,1)�ɈႢ�Ȃ��B

! nb(j) j�Ԗڂ̊ǂɐڑ�����Ă���          �ߓ_�ԍ�
! ne(j) j�Ԗڂ̊ǂɐڑ�����Ă�������Е��̐ߓ_�ԍ�

! ind(n) n�Ԗڂ̐ߓ_��1�F�v�Z�ߓ_�A0�F�����n(const head)�Ȃ̂Ōv�Z���Ȃ��Ă���

! ram(j) j�Ԗڂ̊ǂ̌W���@�w�b�h��=ram(j)*Q**2 ��h=ram*Q**2

!   H(n) n�Ԗڂ̐ߓ_����ђ����r��head, Reservoir(lake) is constant(Boundary condition)


ntn=0; nb=0; ne=0; ram=0.0; nmax=0

 open(10,file='tubes.txt')  ! tube data �Ǐ��i�ߓ_���ł͂Ȃ��j
1 continue
 read(10,*,end=10) n,nb(n),ne(n),ram(n) ! node number of both ends of the tube(n), tube constant ram(n)
 nmax=max(nmax,nb(n))
 nmax=max(nmax,ne(n))
 goto 1

10 close(10)
   nl=n     ! �ǂ̖{��   total tube number
  non=nmax  ! �ߓ_�̑��� total node number

  write(*,*) 'total node number : ',non,' total tube number : ',nl


 do m=1,non ! �ߓ_�ő|�߂���B�ߓ_�����͊��m���H

 is=1
 do n=1,nl ! �ǔԍ��ł̑|��

 if(nb(n).eq.m)then
 is=is+1
 ntn(m,is)=ne(n)   ! �ߓ_m�ƂȂ����Ă���ߓ_�ԍ�
 ntl(m,is)=n       ! ���̊ǔԍ�
! write(*,*) m,ne(n) ; pause
 endif

 if(ne(n).eq.m)then
 is=is+1
 ntn(m,is)=nb(n)   ! �ߓ_m�ƂȂ����Ă���ߓ_�ԍ�
 ntl(m,is)=n       ! ���̊ǔԍ�
! write(*,*) m,nb(n) ; pause
 endif

 enddo ! �ǂł̑|��

 ntn(m,1)=is-1 ! �ߓ_m��ntn(m,1)�̐ߓ_�ƂȂ����Ă���B
 ntl(m,1)=is-1 ! ���̊ǂȂ̂��Bntn(m,1)�Ɠ���

! write(*,*) m,ntn(m,1) ; pause

 enddo ! �ߓ_�ł̑|�� all node

 open(10,file='tnumber.dat')
 do n=1,non
!write(*, '(10i5)') n,ntn(n,1),(ntn(n,j),j=2,ntn(n,1)+1) ; pause
 write(10,'(10i5)') n,ntn(n,1),(ntn(n,j),j=2,ntn(n,1)+1)
 write(10,'(10i5)') n,ntl(n,1),(ntl(n,j),j=2,ntl(n,1)+1)
 enddo
 close(10)

 ind=1 ! �v�Z���ׂ��ߓ_  1 for all nodes but lake node should be 0 because lake water level is constant.

 open(10,file='resor.dat') ! reservoir(lake) data
 read(10,*) nor ! total lake number �����n�����i�ߓ_�̂����̂ǂꂩ�jnor��non�̓���
 sum=0.0
 do n=1,nor
 read(10,*) k,H(k); ind(k)=0 ! �ߓ_k�Ԃ̐���H(k)�̃f�[�^�Aind=0�͗^���邩��v�Z���Ȃ�
 sum=sum+H(k)
 enddo
 ave=sum/nor ! �����r�̕��ϐ��� ! average lake water level -> initial head of all nodes except lake
 close(10)

 do m=1,non
 if(ind(m).eq.1) H(m)=ave ! ���ϐ��ʂ�S�ߓ_(�����r�łȂ�)�̊��蓖�Ă�B�����l initial head of each node except lake
 enddo



! caluculation start



! �����܂ŁA�ߓ_�ԍ�m�ɂ��� ind(m)=1�̐���H(m)�ɂ��Ĕ��W�v�Z���s���Bind(m)=0�ł͒����ʂł�������Ƃ��ė^����
! ���̐ߓ_m�ɂ�ntn(m,1)�̕ʂ̐ߓ_���ڑ�����Ă���A���̔ԍ���ntn(m,2)...ntn(m,ntn(m,1)+1)�Ɋi�[����Ă���
! ���̐ߓ_m�ɂ�ntl(m,1)�{�̊ǂ��ڑ�����Ă���A���̔ԍ���ntn(m,2)..ntn(m,ntn(m,1)+1)�Ɋi�[����Ă���Bntn(m,1)=ntl(m,1)
! �S�ߓ_����non�ł���B
! �ߓ_m�Ɛڑ�����Ă���ʂ̐ߓ_�i�ԍ�ntn(m,n)�j�Ƃ̊Ԃ̊ǔԍ���ntl(m,n)�ł���B


 do it=1,itmax ! �J��Ԃ��v�Z����Biteration

 echeck=0.0      ! �J��Ԃ��v�Z�̎����`�F�b�J�[ convergence checker

 do n=1,non ! �S�ߓ_��|�� all nodes

 sum=0.0    ! �a���[���ɂ����� ��Q initial value

 if(ind(n).eq.1) then ! 1�̏ꍇ�������ׂ�

 do m=1,ntn(n,1) ! nodes connected to node(n)

 i=ntn(n,1+m) ! ����̐ߓ_�ԍ� the node number
 j=ntl(n,1+m) ! ���̊ǔԍ�     the tube number
 dh=H(i)-H(n) ! ���̒l�����̂Ƃ�����{ head difference between i,n
 Q=sqrt(abs(dh)/ram(j)) ! absolute value of flow quantitity no direction

! if conceptual Q direction is from n->i. if i->n ( H(i)>H(n);dh>0 ), then -Q

  if(dh.gt.0.0) Q=-Q ! Q�͏o�Ă�����(n->i)�����Ƃ��߂��̂�H(n)>H(i)����{�B���̋t

 sum=sum+Q    ! �����Ă��� suming 

 enddo ! �ߓ_n�Ɋ֌W���Ă���ʂ̐ߓ_�E�ǘH�ő|��

 echeck=amax1(abs(sum),echeck) ! sum�̓[���ł����ė~�������� ���̍ő�l�𒲂ׂ�

! if(sum>0) too much going out, H(n) goes down
! if(sum<0) too much coming in, H(n) goew up

 H(n)=H(n)-alfa*sum ! sum�̒l�ɉ����ĕύX����Bsum>0�Ȃ猸�炵�Asum<0�Ȃ瑝�₷�B

 endif

 enddo ! �S�ߓ_ for all nodes

if(echeck.lt.emin) then ! �����A�v�Z�I���A�ŏI�����@���ʂƗ��ʋL�^

! just output the results

open(10,file='output.dat')

do m=1,nl ! �ǂ��Ƃɏo��
 h1=H(nb(m)); h2=H(ne(m)); dh=h1-h2
 Q=sqrt(abs(dh)/ram(m))  
 if(dh.lt.0.0) Q=-Q      
 write(10,'(3i5,2f8.2,f9.4)') m, nb(m), ne(m), h1, h2, Q
enddo

close(10)

stop

endif

enddo ! �����J��Ԃ��v�Z iteration

print *, 'no convergence'

open(10,file='outpute.dat')

do m=1,nl ! �ǂ��Ƃɏo��
 h1=H(nb(m)); h2=H(ne(m)); dh=h1-h2
 Q=sqrt(abs(dh)/ram(m))  
 if(dh.lt.0.0) Q=-Q      
 write(10,'(3i5,2f8.2,f9.4)') m, nb(m), ne(m), h1, h2, Q
enddo

close(10)

stop
end
