	.file	"naive.c"
	.intel_syntax noprefix
	.text
	.p2align 4
	.globl	grid_dimension
	.type	grid_dimension, @function
grid_dimension:
.LFB24:
	.cfi_startproc
	push	r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	push	rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	push	rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	mov	rbx, rdx
	lea	rbp, 12[rbx]
	lea	r12, 8[rbx]
	sub	rsp, 16
	.cfi_def_cfa_offset 48
	mov	eax, DWORD PTR [rdi]
	mov	edx, DWORD PTR 4[rdi]
	mov	ecx, eax
	cmp	eax, edx
	jg	.L2
	mov	rcx, r12
	mov	r12, rbp
	mov	rbp, rcx
	mov	ecx, edx
.L2:
	imul	eax, edx
	pxor	xmm0, xmm0
	pxor	xmm2, xmm2
	pxor	xmm1, xmm1
	cvtsi2ss	xmm1, ecx
	cdq
	idiv	esi
	cvtsi2ss	xmm0, eax
	ucomiss	xmm2, xmm0
	ja	.L14
	sqrtss	xmm0, xmm0
.L5:
	divss	xmm1, xmm0
	cvtss2sd	xmm1, xmm1
	addsd	xmm1, QWORD PTR .LC1[rip]
	cvttsd2si	ecx, xmm1
	jmp	.L17
	.p2align 4,,10
	.p2align 3
.L19:
	sub	ecx, 1
.L17:
	mov	eax, esi
	cdq
	idiv	ecx
	test	edx, edx
	jne	.L19
	mov	eax, esi
	mov	DWORD PTR [r12], ecx
	cdq
	idiv	ecx
	mov	DWORD PTR 0[rbp], eax
	mov	QWORD PTR [rbx], rdi
	add	rsp, 16
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	pop	rbx
	.cfi_def_cfa_offset 24
	pop	rbp
	.cfi_def_cfa_offset 16
	pop	r12
	.cfi_def_cfa_offset 8
	ret
.L14:
	.cfi_restore_state
	mov	DWORD PTR 12[rsp], esi
	mov	QWORD PTR [rsp], rdi
	movss	DWORD PTR 8[rsp], xmm1
	call	sqrtf@PLT
	mov	esi, DWORD PTR 12[rsp]
	mov	rdi, QWORD PTR [rsp]
	movss	xmm1, DWORD PTR 8[rsp]
	jmp	.L5
	.cfi_endproc
.LFE24:
	.size	grid_dimension, .-grid_dimension
	.p2align 4
	.globl	get_frame
	.type	get_frame, @function
get_frame:
.LFB25:
	.cfi_startproc
	push	r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	mov	r9d, edx
	mov	rax, rdi
	mov	r11d, esi
	push	r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	push	r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	mov	r12d, ecx
	push	rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	xor	ebp, ebp
	push	rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	mov	rdx, QWORD PTR [rdi]
	mov	ebx, DWORD PTR 8[rax]
	mov	ecx, DWORD PTR 12[rax]
	mov	QWORD PTR [r8], rdi
	mov	edi, DWORD PTR [rdx]
	mov	esi, DWORD PTR 4[rdx]
	mov	eax, edi
	cdq
	idiv	ebx
	mov	r10d, eax
	mov	eax, esi
	cdq
	idiv	ecx
	mov	edx, ebx
	imul	ebx, r10d
	sub	edx, edi
	lea	edx, -1[rdx+rbx]
	mov	ebx, ecx
	sub	ebx, esi
	imul	ecx, eax
	cmp	edx, r12d
	setl	bpl
	add	ebp, r10d
	lea	r14d, -1[rbx+rcx]
	xor	ebx, ebx
	mov	ecx, r12d
	cmp	r14d, r9d
	setl	bl
	sub	ecx, edx
	xor	r13d, r13d
	add	ebx, eax
	sub	ecx, 1
	mov	edx, ecx
	cmovs	edx, r13d
	imul	r10d, r12d
	lea	ecx, [rdx+r10]
	mov	edx, r9d
	mov	r10d, r11d
	sub	edx, r14d
	sub	edx, 1
	cmovs	edx, r13d
	imul	eax, r9d
	mov	r9d, r11d
	add	edx, eax
	cmp	ecx, r11d
	mov	eax, edi
	cmovle	r9d, ecx
	cmp	edx, r11d
	cmovle	r10d, edx
	sub	eax, ebp
	sub	eax, ecx
	cmp	eax, r11d
	cmovg	eax, r11d
	sub	esi, ebx
	sub	esi, edx
	add	eax, r9d
	cmp	esi, r11d
	cmovg	esi, r11d
	sub	edx, r10d
	sub	ecx, r9d
	add	eax, ebp
	imul	edx, edi
	mov	DWORD PTR 8[r8], eax
	add	esi, r10d
	add	esi, ebx
	pop	rbx
	.cfi_def_cfa_offset 40
	pop	rbp
	.cfi_def_cfa_offset 32
	add	edx, ecx
	pop	r12
	.cfi_def_cfa_offset 24
	pop	r13
	.cfi_def_cfa_offset 16
	mov	DWORD PTR 12[r8], esi
	mov	DWORD PTR 16[r8], edx
	pop	r14
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE25:
	.size	get_frame, .-get_frame
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC2:
	.string	"[THREAD %d] coordinates: (%d, %d)\n            dimension: (%d, %d)\n            p_start: %d\n"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"threads number: %d\n"
.LC5:
	.string	"allocated memory: %f Mb.\n"
	.section	.rodata.str1.8
	.align 8
.LC6:
	.string	"w_threads: %d, h_threads: %d.\n\n"
	.text
	.p2align 4
	.type	main._omp_fn.0, @function
main._omp_fn.0:
.LFB27:
	.cfi_startproc
	push	r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	push	r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	mov	rbp, rdi
	push	rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	sub	rsp, 72
	.cfi_def_cfa_offset 128
	mov	r13d, DWORD PTR 24[rdi]
	mov	rax, QWORD PTR fs:40
	mov	QWORD PTR 56[rsp], rax
	xor	eax, eax
	call	omp_get_thread_num@PLT
	mov	r12d, eax
	test	eax, eax
	je	.L23
.L33:
	call	GOMP_barrier@PLT
	mov	rsi, QWORD PTR 8[rbp]
	mov	eax, r12d
	lea	rdi, .gomp_critical_user_show[rip]
	cdq
	idiv	DWORD PTR 8[rsi]
	mov	r14d, eax
	mov	r15d, edx
	call	GOMP_critical_name_start@PLT
	mov	rdi, QWORD PTR 8[rbp]
	mov	ecx, r15d
	mov	edx, r14d
	lea	r8, 32[rsp]
	mov	esi, r13d
	call	get_frame
	sub	rsp, 8
	.cfi_def_cfa_offset 136
	mov	edx, r14d
	mov	ecx, r15d
	mov	eax, DWORD PTR 56[rsp]
	mov	esi, r12d
	lea	rdi, .LC2[rip]
	push	rax
	.cfi_def_cfa_offset 144
	mov	r9d, DWORD PTR 60[rsp]
	xor	eax, eax
	mov	r8d, DWORD PTR 56[rsp]
	call	printf@PLT
	lea	rdi, .gomp_critical_user_show[rip]
	call	GOMP_critical_name_end@PLT
	pop	rax
	.cfi_def_cfa_offset 136
	pop	rdx
	.cfi_def_cfa_offset 128
	mov	rax, QWORD PTR 56[rsp]
	sub	rax, QWORD PTR fs:40
	jne	.L44
	add	rsp, 72
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	pop	rbx
	.cfi_def_cfa_offset 48
	pop	rbp
	.cfi_def_cfa_offset 40
	pop	r12
	.cfi_def_cfa_offset 32
	pop	r13
	.cfi_def_cfa_offset 24
	pop	r14
	.cfi_def_cfa_offset 16
	pop	r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L23:
	.cfi_restore_state
	call	omp_get_num_threads@PLT
	mov	ebx, DWORD PTR 20[rbp]
	mov	r14d, DWORD PTR 16[rbp]
	lea	rdi, .LC3[rip]
	mov	DWORD PTR 28[rbp], eax
	mov	esi, eax
	xor	eax, eax
	imul	ebx, r14d
	call	printf@PLT
	movsx	rax, ebx
	add	rax, rax
	js	.L24
	pxor	xmm0, xmm0
	cvtsi2ss	xmm0, rax
.L25:
	mulss	xmm0, DWORD PTR .LC4[rip]
	lea	rdi, .LC5[rip]
	mov	eax, 1
	cvtss2sd	xmm0, xmm0
	call	printf@PLT
	mov	r15, QWORD PTR 0[rbp]
	mov	r14, QWORD PTR 8[rbp]
	mov	ebx, DWORD PTR 28[rbp]
	mov	eax, DWORD PTR [r15]
	mov	edx, DWORD PTR 4[r15]
	cmp	eax, edx
	jle	.L26
	lea	rdi, 8[r14]
	lea	rsi, 12[r14]
	mov	ecx, eax
.L27:
	imul	eax, edx
	pxor	xmm0, xmm0
	pxor	xmm2, xmm2
	pxor	xmm1, xmm1
	cvtsi2ss	xmm1, ecx
	cdq
	idiv	ebx
	cvtsi2ss	xmm0, eax
	ucomiss	xmm2, xmm0
	ja	.L41
	sqrtss	xmm0, xmm0
	mov	r8, r14
.L30:
	divss	xmm1, xmm0
	cvtss2sd	xmm1, xmm1
	addsd	xmm1, QWORD PTR .LC1[rip]
	cvttsd2si	ecx, xmm1
	jmp	.L43
	.p2align 4,,10
	.p2align 3
.L45:
	sub	ecx, 1
.L43:
	mov	eax, ebx
	cdq
	idiv	ecx
	test	edx, edx
	jne	.L45
	mov	eax, ebx
	mov	DWORD PTR [rdi], ecx
	lea	rdi, .LC6[rip]
	cdq
	idiv	ecx
	mov	DWORD PTR [rsi], eax
	mov	edx, DWORD PTR 12[r8]
	xor	eax, eax
	mov	QWORD PTR [r14], r15
	mov	esi, DWORD PTR 8[r8]
	call	printf@PLT
	jmp	.L33
	.p2align 4,,10
	.p2align 3
.L24:
	shr	rax
	pxor	xmm0, xmm0
	cvtsi2ss	xmm0, rax
	addss	xmm0, xmm0
	jmp	.L25
	.p2align 4,,10
	.p2align 3
.L26:
	lea	rdi, 12[r14]
	lea	rsi, 8[r14]
	mov	ecx, edx
	jmp	.L27
.L41:
	mov	QWORD PTR 16[rsp], rsi
	mov	QWORD PTR 8[rsp], rdi
	movss	DWORD PTR 28[rsp], xmm1
	call	sqrtf@PLT
	mov	r8, QWORD PTR 8[rbp]
	movss	xmm1, DWORD PTR 28[rsp]
	mov	rsi, QWORD PTR 16[rsp]
	mov	rdi, QWORD PTR 8[rsp]
	jmp	.L30
.L44:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE27:
	.size	main._omp_fn.0, .-main._omp_fn.0
	.section	.rodata.str1.8
	.align 8
.LC7:
	.string	"USAGE:\n%s <width> <height> <kernel_size>\n"
	.align 8
.LC8:
	.string	"Bad allocation: not enough memory.\n"
	.section	.rodata.str1.1
.LC9:
	.string	"PATH"
.LC11:
	.string	"Data initialized by master."
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB26:
	.cfi_startproc
	push	r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	push	r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	mov	rbx, rsi
	sub	rsp, 120
	.cfi_def_cfa_offset 176
	mov	DWORD PTR 12[rsp], edi
	mov	rcx, QWORD PTR fs:40
	mov	QWORD PTR 104[rsp], rcx
	xor	ecx, ecx
	cmp	edi, 4
	je	.L47
	mov	rdx, QWORD PTR [rsi]
	mov	rdi, QWORD PTR stderr[rip]
	xor	eax, eax
	lea	rsi, .LC7[rip]
	call	fprintf@PLT
	mov	eax, 1
.L46:
	mov	rcx, QWORD PTR 104[rsp]
	sub	rcx, QWORD PTR fs:40
	jne	.L57
	add	rsp, 120
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	pop	rbx
	.cfi_def_cfa_offset 48
	pop	rbp
	.cfi_def_cfa_offset 40
	pop	r12
	.cfi_def_cfa_offset 32
	pop	r13
	.cfi_def_cfa_offset 24
	pop	r14
	.cfi_def_cfa_offset 16
	pop	r15
	.cfi_def_cfa_offset 8
	ret
.L47:
	.cfi_restore_state
	mov	rdi, QWORD PTR 8[rsi]
	mov	edx, 10
	xor	esi, esi
	call	strtol@PLT
	mov	rdi, QWORD PTR 16[rbx]
	mov	edx, 10
	xor	esi, esi
	mov	r13, rax
	call	strtol@PLT
	mov	rdi, QWORD PTR 24[rbx]
	mov	edx, 10
	xor	esi, esi
	mov	r12, rax
	call	strtol@PLT
	movd	xmm3, r12d
	imul	r12d, r13d
	movd	xmm2, r13d
	punpckldq	xmm2, xmm3
	mov	esi, 2
	mov	r14, rax
	movq	QWORD PTR [rsp], xmm2
	movsx	rdi, r12d
	movq	QWORD PTR 24[rsp], xmm2
	call	calloc@PLT
	mov	r13, rax
	test	rax, rax
	je	.L58
	lea	rdi, .LC9[rip]
	lea	rbp, 80[rsp]
	call	getenv@PLT
	lea	r15, 48[rsp]
	mov	rdi, rax
	call	strlen@PLT
	mov	rbx, rax
	call	getpid@PLT
	mov	rsi, rbp
	movsx	rdi, eax
	lea	rax, 12[rsp]
	imul	rdi, rbx
	xor	rdi, rax
	call	srand48_r@PLT
	test	r12d, r12d
	jle	.L50
	lea	eax, -1[r12]
	mov	rbx, r13
	lea	r15, 48[rsp]
	add	rax, rax
	lea	r12, 2[r13+rax]
	.p2align 4,,10
	.p2align 3
.L51:
	mov	rsi, r15
	mov	rdi, rbp
	add	rbx, 2
	call	drand48_r@PLT
	movsd	xmm0, QWORD PTR .LC10[rip]
	mulsd	xmm0, QWORD PTR 48[rsp]
	addsd	xmm0, QWORD PTR .LC1[rip]
	cvttsd2si	eax, xmm0
	mov	WORD PTR -2[rbx], ax
	cmp	rbx, r12
	jne	.L51
.L50:
	lea	rdi, .LC11[rip]
	call	puts@PLT
	lea	rax, 24[rsp]
	xor	ecx, ecx
	xor	edx, edx
	movq	xmm0, rax
	lea	rax, 32[rsp]
	movd	xmm1, r14d
	mov	rsi, r15
	movq	xmm4, rax
	lea	rdi, main._omp_fn.0[rip]
	punpcklqdq	xmm0, xmm4
	movaps	XMMWORD PTR 48[rsp], xmm0
	movq	xmm0, QWORD PTR [rsp]
	punpcklqdq	xmm0, xmm1
	movaps	XMMWORD PTR 64[rsp], xmm0
	call	GOMP_parallel@PLT
	mov	rdi, r13
	call	free@PLT
	xor	eax, eax
	jmp	.L46
.L57:
	call	__stack_chk_fail@PLT
.L58:
	mov	rcx, QWORD PTR stderr[rip]
	mov	edx, 35
	mov	esi, 1
	lea	rdi, .LC8[rip]
	call	fwrite@PLT
	mov	eax, 1
	jmp	.L46
	.cfi_endproc
.LFE26:
	.size	main, .-main
	.comm	.gomp_critical_user_show,8,8
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC1:
	.long	0
	.long	1071644672
	.section	.rodata.cst4,"aM",@progbits,4
	.align 4
.LC4:
	.long	897581056
	.section	.rodata.cst8
	.align 8
.LC10:
	.long	0
	.long	1089470432
	.ident	"GCC: (GNU) 10.2.0"
	.section	.note.GNU-stack,"",@progbits
