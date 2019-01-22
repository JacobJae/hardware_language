/* 
--------------------------------------------------------
Title : CPU
Name: Sangheon_Jae ID: 213166392
--------------------------------------------------------
*/

module yMux1(z, a, b, c);

output z;
input a, b, c;
wire notC, upper, lower;

not my_not(notC, c);
and upperAnd(upper, a, notC);
and lowerAnd(lower, c, b);
or  my_or(z, upper, lower);

endmodule
////////////////////////////////////////////////////////

module yMux(z, a, b, c);

parameter SIZE = 2;
output [SIZE-1:0] z;
input [SIZE-1:0] a, b;
input c;

yMux1 mine[SIZE-1:0](z, a, b, c);

endmodule
////////////////////////////////////////////////////////

module yMux4to1(z, a0,a1,a2,a3, c);

parameter SIZE = 32;
output [SIZE-1:0] z;
input [SIZE-1:0] a0, a1, a2, a3;
input [1:0] c;
wire [SIZE-1:0] zLo, zHi;

yMux #(SIZE) lo(zLo, a0, a1, c[0]);
yMux #(SIZE) hi(zHi, a2, a3, c[0]);
yMux #(SIZE) final(z, zLo, zHi, c[1]);

endmodule
////////////////////////////////////////////////////////

module yAdder1(z, cout, a, b, cin);

output z, cout;
input a, b, cin;

xor left_xor(tmp, a, b);
xor right_xor(z, cin, tmp);
and left_and(outL, a, b);
and right_and(outR, tmp, cin);
or my_or(cout, outR, outL);

endmodule
////////////////////////////////////////////////////////

module yAdder(z, cout, a, b, cin);

output [31:0] z;
output cout;
input [31:0] a, b;
input cin;
wire[31:0] in, out;

yAdder1 mine[31:0](z, out, a, b, in);
assign in[0] = cin;
assign in[31:1] = out[30:0];
assign cout = out[31];

endmodule
////////////////////////////////////////////////////////

module yArith(z, cout, a, b, ctrl);

// add if ctrl=0, subtract if ctrl=1
output [31:0] z;
output cout;
input [31:0] a, b;
input ctrl;
wire[31:0] notB, tmp;
wire cin;

not my_not[31:0](notB, b);
assign tmp = (ctrl) ? notB : b;
yAdder mine(z, out, a, tmp, ctrl);

endmodule
////////////////////////////////////////////////////////

module yAlu(z, zero, a, b, op);

input [31:0] a, b;
input [2:0] op;
output [31:0] z;
output zero;
wire ex1, ex2, z1;
wire [31:0] i0, i1, i2, i3, slt;
wire [15:0] z16;
wire [7:0] z8;
wire [3:0] z4;
wire [1:0] z2;
assign slt[31:1] = 0;
assign ex = 0; // not supported

yArith my_subs(i3, ex2, a, b, 1'b1);
assign slt[0] = (a[31] != b[31]) ? a[31] : i3[31];

and my_and[31:0](i0, a, b);
or  my_or [31:0](i1, a, b);
yArith my_arith(i2, ex1, a, b, op[2]);
yMux4to1 my_mux(z, i0, i1, i2, slt, op[1:0]);
or or16[15:0] (z16, z[15:0], z[31:16]);
or or8[7:0] (z8, z16[7:0], z16[15:8]); 
or or4[3:0] (z4, z8[3:0], z8[7:4]);
or or2[1:0] (z2, z4[1:0], z4[3:2]);
or or1		(z1, z2[1], z2[0]);
assign zero = ~z1;

endmodule
////////////////////////////////////////////////////////

module yIF(ins, PCp4, PCin, clk);

	output 	[31:0] 	ins, PCp4;
	input	[31:0] 	PCin;
	input 		clk;
	reg	[31:0]	PC;
	wire		zero;

	yAlu my_Alu (PCp4, zero, PC, 4, 3'b010);
	mem my_mem (ins, PC, 0, clk, 1'b1, 1'b0);

	always@(posedge clk)
	begin
		PC <= PCin;
	end
endmodule
////////////////////////////////////////////////////////

module yID(rd1, rd2, imm, jTarget, ins, wd, RegDst, RegWrite, clk);

	output 	[31:0] 	rd1, rd2, imm;
	output	[25:0] 	jTarget;
	input 	[31:0]	ins, wd;
	input 		RegDst, RegWrite, clk;
	wire	[4:0]	wn, rn1, rn2;

	assign rn1 = ins[25:21];
	assign rn2 = ins[20:16];
	yMux #(5) my_mux1 (wn, rn2, ins[15:11], RegDst); 
	rf myRF(rd1, rd2, rn1, rn2, wn, wd, clk, RegWrite);
	assign jTarget = ins[25:0];
	assign imm[15:0] = ins[15:0];
	yMux #(16) my_mux2 (imm[31:16], 16'b0000000000000000, 16'b1111111111111111, ins[15]); 
endmodule
////////////////////////////////////////////////////////

module yEX(z, zero, rd1, rd2, imm, op, ALUSrc);

output 	[31:0] 	z;
output 		zero;
input	[31:0]	rd1, rd2, imm;
input 	[2:0]	op;
input 		ALUSrc; 
wire	[31:0]	b;

	yMux #(32) my_mux1 (b, rd2, imm, ALUSrc); 
	yAlu my_Alu (z, zero, rd1, b, op);
endmodule
////////////////////////////////////////////////////////

module yDM(memOut, exeOut, rd2, clk, MemRead, MemWrite);

output [31:0] memOut;
input [31:0] exeOut, rd2;
input clk, MemRead, MemWrite;

// instantiate the circuit (only one line)
mem my_mem (memOut, exeOut, rd2, clk, MemRead, MemWrite);

endmodule
//////////////////////////////////////////////////////////

module yWB(wb, exeOut, memOut, Mem2Reg);

output [31:0] wb;
input [31:0] exeOut, memOut;
input Mem2Reg;

// instantiate the circuit (only one line)
yMux #(32) my_mux (wb, exeOut, memOut, Mem2Reg);

endmodule
////////////////////////////////////////////////////////

module yPC(PCin, PCp4,INT,entryPoint,imm,jTarget,zero,branch,jump);

output [31:0] PCin;
input [31:0] PCp4, entryPoint, imm;
input [25:0] jTarget;
input zero, branch, jump, INT;
wire [31:0] immX4, jimm, bTarget, choiceA, choiceB;
wire doBranch, zf;

assign immX4[31:2] = imm[29:0];
assign immX4[1:0] = 2'b00;
assign jimm[31:28] = PCp4[31:28];
assign jimm[27:2] = jTarget;
assign jimm[1:0] = 2'b00;

yAlu myALU(bTarget, zf, PCp4, immX4, 3'b010);
and (doBranch, branch, zero);
yMux #(32) mux1(choiceA, PCp4, bTarget, doBranch);
yMux #(32) mux2(choiceB, choiceA, jimm, jump);
yMux #(32) mux3(PCin, choiceB, entryPoint, INT);

endmodule
////////////////////////////////////////////////////////

module yC1(rtype, lw, sw, jump, branch, opCode);

output rtype, lw, sw, jump, branch;
input [5:0] opCode;
wire not4, not3, not2, not1;

not (not4, opCode[4]);
not (not3, opCode[3]);
not (not2, opCode[2]);
not (not1, opCode[1]);

// generate lw, opCode is 100011
and (lw, opCode[5], not4, not3, not2, opCode[1], opCode[0]);

// generate sw, opCode is 101011
and (sw, opCode[5], not4, opCode[3], not2, opCode[1], opCode[0]);

// generate branch, opCode is 000100
nor (branch, opCode[5], opCode[4], opCode[3], not2, opCode[1], opCode[0]);

// generate jump, opCode is 000010
nor (jump, opCode[5], opCode[4], opCode[3], opCode[2], not1, opCode[0]);

// generate R-type, opCode is 000000
nor (rtype, opCode[5], opCode[4], opCode[3], opCode[2], opCode[1], opCode[0]);

endmodule
////////////////////////////////////////////////////////

module yC2(RegDst, ALUSrc, RegWrite, Mem2Reg, MemRead, MemWrite, rtype, lw, sw, branch);

output RegDst, ALUSrc, RegWrite, Mem2Reg, MemRead, MemWrite;
input rtype, lw, sw, branch;

assign RegDst = rtype;
assign Mem2Reg = lw;
assign MemRead = lw;
assign MemWrite = sw;

nor (ALUSrc, rtype, branch);
nor (RegWrite, sw, branch);
endmodule
////////////////////////////////////////////////////////

module yC3(ALUop, rtype, branch);

output [1:0] ALUop;
input rtype, branch;

assign ALUop[1] = rtype;
assign ALUop[0] = branch;

endmodule
////////////////////////////////////////////////////////

module yC4(op, ALUop, fnCode);

output [2:0] op;
input [5:0] fnCode;
input [1:0] ALUop;

// instantiate and connect
or (or1, fnCode[0], fnCode[3]);
and (and1, ALUop[1], fnCode[1]);
or (op[2], and1, ALUop[0]);
nand (op[1], fnCode[2], ALUop[1]);
and (op[0], or1, ALUop[1]);

endmodule 
////////////////////////////////////////////////////////

module yChip(ins, rd2, wb, entryPoint, INT, clk);

output [31:0] ins, rd2, wb;
input [31:0] entryPoint;
input INT, clk;
wire [31:0] PCin, wd, wb, memOut, rd1, rd2, imm, ins, PCp4, z;
wire [25:0] jTarget;
wire [5:0] opCode, fnCode;
wire [2:0] op;
wire [1:0] ALUop;
wire zero, rtype, lw, sw, jump, branch, RegDst, RegWrite, ALUSrc, MemRead, MemWrite, Mem2Reg;

yIF myIF(ins, PCp4, PCin, clk);
yID myID(rd1, rd2, imm, jTarget, ins, wd, RegDst, RegWrite, clk);
yEX myEx(z, zero, rd1, rd2, imm, op, ALUSrc);
yDM myDM(memOut, z, rd2, clk, MemRead, MemWrite);
yWB myWB(wb, z, memOut, Mem2Reg);
assign wd = wb; 
yPC myPC(PCin, PCp4, INT, entryPoint, imm, jTarget, zero, branch, jump);

assign opCode = ins[31:26];
yC1 myC1(rtype, lw, sw, jump, branch, opCode);
yC2 myC2(RegDst, ALUSrc, RegWrite, Mem2Reg, MemRead, MemWrite, rtype, lw, sw, branch); 

assign fnCode = ins[5:0];
yC3 myC3(ALUop, rtype, branch);
yC4 myC4(op, ALUop, fnCode);

endmodule 
////////////////////////////////////////////////////////


////////////////////////////////////////////////////////