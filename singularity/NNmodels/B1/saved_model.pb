ът
пД
,
Abs
x"T
y"T"
Ttype:

2	
:
Add
x"T
y"T
z"T"
Ttype:
2	
E
AssignAddVariableOp
resource
value"dtype"
dtypetype
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
8
Const
output"dtype"
valuetensor"
dtypetype
B
GreaterEqual
x"T
y"T
z
"
Ttype:
2	
.
Identity

input"T
output"T"	
Ttype
,
Log
x"T
y"T"
Ttype:

2
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
8
Maximum
x"T
y"T
z"T"
Ttype:

2	

Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(
8
Minimum
x"T
y"T
z"T"
Ttype:

2	
=
Mul
x"T
y"T
z"T"
Ttype:
2	

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
>
RealDiv
x"T
y"T
z"T"
Ttype:
2	
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
.
Rsqrt
x"T
y"T"
Ttype:

2
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
H
ShardedFilename
basename	
shard

num_shards
filename
O
Size

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
@
Softplus
features"T
activations"T"
Ttype:
2
-
Sqrt
x"T
y"T"
Ttype:

2
1
Square
x"T
y"T"
Ttype:

2	
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
:
Sub
x"T
y"T
z"T"
Ttype:
2	

Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
-
Tanh
x"T
y"T"
Ttype:

2
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"serve*1.14.02v1.14.0-rc1-22-gaf24dc9Ч
r
dense_1_inputPlaceholder*
dtype0*(
_output_shapes
:џџџџџџџџџ*
shape:џџџџџџџџџ
m
dense_1/random_uniform/shapeConst*
valueB"      *
dtype0*
_output_shapes
:
_
dense_1/random_uniform/minConst*
valueB
 *ѓ5Н*
dtype0*
_output_shapes
: 
_
dense_1/random_uniform/maxConst*
valueB
 *ѓ5=*
dtype0*
_output_shapes
: 
Њ
$dense_1/random_uniform/RandomUniformRandomUniformdense_1/random_uniform/shape*
T0*
dtype0* 
_output_shapes
:
*
seed2јъч*
seedБџх)
z
dense_1/random_uniform/subSubdense_1/random_uniform/maxdense_1/random_uniform/min*
T0*
_output_shapes
: 

dense_1/random_uniform/mulMul$dense_1/random_uniform/RandomUniformdense_1/random_uniform/sub* 
_output_shapes
:
*
T0

dense_1/random_uniformAdddense_1/random_uniform/muldense_1/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_1/kernelVarHandleOp*
shape:
*
dtype0*
_output_shapes
: *
shared_namedense_1/kernel*!
_class
loc:@dense_1/kernel*
	container 
m
/dense_1/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/kernel*
_output_shapes
: 

dense_1/kernel/AssignAssignVariableOpdense_1/kerneldense_1/random_uniform*!
_class
loc:@dense_1/kernel*
dtype0

"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*!
_class
loc:@dense_1/kernel*
dtype0* 
_output_shapes
:

\
dense_1/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
Ѓ
dense_1/biasVarHandleOp*
shared_namedense_1/bias*
_class
loc:@dense_1/bias*
	container *
shape:*
dtype0*
_output_shapes
: 
i
-dense_1/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/bias*
_output_shapes
: 
r
dense_1/bias/AssignAssignVariableOpdense_1/biasdense_1/Const*
_class
loc:@dense_1/bias*
dtype0

 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
_class
loc:@dense_1/bias*
dtype0*
_output_shapes	
:
n
dense_1/MatMul/ReadVariableOpReadVariableOpdense_1/kernel*
dtype0* 
_output_shapes
:


dense_1/MatMulMatMuldense_1_inputdense_1/MatMul/ReadVariableOp*
T0*(
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( 
h
dense_1/BiasAdd/ReadVariableOpReadVariableOpdense_1/bias*
dtype0*
_output_shapes	
:

dense_1/BiasAddBiasAdddense_1/MatMuldense_1/BiasAdd/ReadVariableOp*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ*
T0
\
keras_learning_phase/inputConst*
value	B
 Z *
dtype0
*
_output_shapes
: 
|
keras_learning_phasePlaceholderWithDefaultkeras_learning_phase/input*
dtype0
*
_output_shapes
: *
shape: 
n
dropout_1/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_1/cond/switch_tIdentitydropout_1/cond/Switch:1*
T0
*
_output_shapes
: 
[
dropout_1/cond/switch_fIdentitydropout_1/cond/Switch*
_output_shapes
: *
T0

Y
dropout_1/cond/pred_idIdentitykeras_learning_phase*
T0
*
_output_shapes
: 
z
dropout_1/cond/dropout/rateConst^dropout_1/cond/switch_t*
valueB
 *ЭЬЬ=*
dtype0*
_output_shapes
: 

dropout_1/cond/dropout/ShapeShape%dropout_1/cond/dropout/Shape/Switch:1*
T0*
out_type0*
_output_shapes
:
С
#dropout_1/cond/dropout/Shape/SwitchSwitchdense_1/BiasAdddropout_1/cond/pred_id*
T0*"
_class
loc:@dense_1/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

)dropout_1/cond/dropout/random_uniform/minConst^dropout_1/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_1/cond/dropout/random_uniform/maxConst^dropout_1/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
С
3dropout_1/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_1/cond/dropout/Shape*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2Нпл*
seedБџх)*
T0
Ї
)dropout_1/cond/dropout/random_uniform/subSub)dropout_1/cond/dropout/random_uniform/max)dropout_1/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
У
)dropout_1/cond/dropout/random_uniform/mulMul3dropout_1/cond/dropout/random_uniform/RandomUniform)dropout_1/cond/dropout/random_uniform/sub*
T0*(
_output_shapes
:џџџџџџџџџ
Е
%dropout_1/cond/dropout/random_uniformAdd)dropout_1/cond/dropout/random_uniform/mul)dropout_1/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_1/cond/dropout/sub/xConst^dropout_1/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_1/cond/dropout/subSubdropout_1/cond/dropout/sub/xdropout_1/cond/dropout/rate*
T0*
_output_shapes
: 

 dropout_1/cond/dropout/truediv/xConst^dropout_1/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

dropout_1/cond/dropout/truedivRealDiv dropout_1/cond/dropout/truediv/xdropout_1/cond/dropout/sub*
T0*
_output_shapes
: 
Њ
#dropout_1/cond/dropout/GreaterEqualGreaterEqual%dropout_1/cond/dropout/random_uniformdropout_1/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_1/cond/dropout/mulMul%dropout_1/cond/dropout/Shape/Switch:1dropout_1/cond/dropout/truediv*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_1/cond/dropout/CastCast#dropout_1/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0

dropout_1/cond/dropout/mul_1Muldropout_1/cond/dropout/muldropout_1/cond/dropout/Cast*
T0*(
_output_shapes
:џџџџџџџџџ
Е
dropout_1/cond/Switch_1Switchdense_1/BiasAdddropout_1/cond/pred_id*
T0*"
_class
loc:@dense_1/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

dropout_1/cond/MergeMergedropout_1/cond/Switch_1dropout_1/cond/dropout/mul_1*
N**
_output_shapes
:џџџџџџџџџ: *
T0
b
activation_1/ReluReludropout_1/cond/Merge*
T0*(
_output_shapes
:џџџџџџџџџ
m
dense_2/random_uniform/shapeConst*
valueB"      *
dtype0*
_output_shapes
:
_
dense_2/random_uniform/minConst*
valueB
 *  Н*
dtype0*
_output_shapes
: 
_
dense_2/random_uniform/maxConst*
valueB
 *  =*
dtype0*
_output_shapes
: 
Њ
$dense_2/random_uniform/RandomUniformRandomUniformdense_2/random_uniform/shape*
T0*
dtype0* 
_output_shapes
:
*
seed2јД*
seedБџх)
z
dense_2/random_uniform/subSubdense_2/random_uniform/maxdense_2/random_uniform/min*
T0*
_output_shapes
: 

dense_2/random_uniform/mulMul$dense_2/random_uniform/RandomUniformdense_2/random_uniform/sub* 
_output_shapes
:
*
T0

dense_2/random_uniformAdddense_2/random_uniform/muldense_2/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_2/kernelVarHandleOp*
	container *
shape:
*
dtype0*
_output_shapes
: *
shared_namedense_2/kernel*!
_class
loc:@dense_2/kernel
m
/dense_2/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_2/kernel*
_output_shapes
: 

dense_2/kernel/AssignAssignVariableOpdense_2/kerneldense_2/random_uniform*
dtype0*!
_class
loc:@dense_2/kernel

"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*!
_class
loc:@dense_2/kernel*
dtype0* 
_output_shapes
:

\
dense_2/ConstConst*
dtype0*
_output_shapes	
:*
valueB*    
Ѓ
dense_2/biasVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_2/bias*
_class
loc:@dense_2/bias*
	container *
shape:
i
-dense_2/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_2/bias*
_output_shapes
: 
r
dense_2/bias/AssignAssignVariableOpdense_2/biasdense_2/Const*
_class
loc:@dense_2/bias*
dtype0

 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_class
loc:@dense_2/bias*
dtype0*
_output_shapes	
:
n
dense_2/MatMul/ReadVariableOpReadVariableOpdense_2/kernel*
dtype0* 
_output_shapes
:

Ѓ
dense_2/MatMulMatMulactivation_1/Reludense_2/MatMul/ReadVariableOp*(
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( *
T0
h
dense_2/BiasAdd/ReadVariableOpReadVariableOpdense_2/bias*
dtype0*
_output_shapes	
:

dense_2/BiasAddBiasAdddense_2/MatMuldense_2/BiasAdd/ReadVariableOp*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ*
T0
n
dropout_2/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_2/cond/switch_tIdentitydropout_2/cond/Switch:1*
T0
*
_output_shapes
: 
[
dropout_2/cond/switch_fIdentitydropout_2/cond/Switch*
T0
*
_output_shapes
: 
Y
dropout_2/cond/pred_idIdentitykeras_learning_phase*
T0
*
_output_shapes
: 
z
dropout_2/cond/dropout/rateConst^dropout_2/cond/switch_t*
valueB
 *ЭЬЬ=*
dtype0*
_output_shapes
: 

dropout_2/cond/dropout/ShapeShape%dropout_2/cond/dropout/Shape/Switch:1*
T0*
out_type0*
_output_shapes
:
С
#dropout_2/cond/dropout/Shape/SwitchSwitchdense_2/BiasAdddropout_2/cond/pred_id*
T0*"
_class
loc:@dense_2/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

)dropout_2/cond/dropout/random_uniform/minConst^dropout_2/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_2/cond/dropout/random_uniform/maxConst^dropout_2/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
С
3dropout_2/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_2/cond/dropout/Shape*
seedБџх)*
T0*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2ѓдЅ
Ї
)dropout_2/cond/dropout/random_uniform/subSub)dropout_2/cond/dropout/random_uniform/max)dropout_2/cond/dropout/random_uniform/min*
_output_shapes
: *
T0
У
)dropout_2/cond/dropout/random_uniform/mulMul3dropout_2/cond/dropout/random_uniform/RandomUniform)dropout_2/cond/dropout/random_uniform/sub*(
_output_shapes
:џџџџџџџџџ*
T0
Е
%dropout_2/cond/dropout/random_uniformAdd)dropout_2/cond/dropout/random_uniform/mul)dropout_2/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_2/cond/dropout/sub/xConst^dropout_2/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_2/cond/dropout/subSubdropout_2/cond/dropout/sub/xdropout_2/cond/dropout/rate*
_output_shapes
: *
T0

 dropout_2/cond/dropout/truediv/xConst^dropout_2/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *  ?

dropout_2/cond/dropout/truedivRealDiv dropout_2/cond/dropout/truediv/xdropout_2/cond/dropout/sub*
_output_shapes
: *
T0
Њ
#dropout_2/cond/dropout/GreaterEqualGreaterEqual%dropout_2/cond/dropout/random_uniformdropout_2/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_2/cond/dropout/mulMul%dropout_2/cond/dropout/Shape/Switch:1dropout_2/cond/dropout/truediv*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_2/cond/dropout/CastCast#dropout_2/cond/dropout/GreaterEqual*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0*

SrcT0


dropout_2/cond/dropout/mul_1Muldropout_2/cond/dropout/muldropout_2/cond/dropout/Cast*
T0*(
_output_shapes
:џџџџџџџџџ
Е
dropout_2/cond/Switch_1Switchdense_2/BiasAdddropout_2/cond/pred_id*
T0*"
_class
loc:@dense_2/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

dropout_2/cond/MergeMergedropout_2/cond/Switch_1dropout_2/cond/dropout/mul_1*
T0*
N**
_output_shapes
:џџџџџџџџџ: 
b
activation_2/ReluReludropout_2/cond/Merge*
T0*(
_output_shapes
:џџџџџџџџџ
m
dense_3/random_uniform/shapeConst*
valueB"      *
dtype0*
_output_shapes
:
_
dense_3/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *ѓЕН
_
dense_3/random_uniform/maxConst*
valueB
 *ѓЕ=*
dtype0*
_output_shapes
: 
Њ
$dense_3/random_uniform/RandomUniformRandomUniformdense_3/random_uniform/shape*
dtype0* 
_output_shapes
:
*
seed2Ьч*
seedБџх)*
T0
z
dense_3/random_uniform/subSubdense_3/random_uniform/maxdense_3/random_uniform/min*
T0*
_output_shapes
: 

dense_3/random_uniform/mulMul$dense_3/random_uniform/RandomUniformdense_3/random_uniform/sub*
T0* 
_output_shapes
:


dense_3/random_uniformAdddense_3/random_uniform/muldense_3/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_3/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_3/kernel*!
_class
loc:@dense_3/kernel*
	container *
shape:

m
/dense_3/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_3/kernel*
_output_shapes
: 

dense_3/kernel/AssignAssignVariableOpdense_3/kerneldense_3/random_uniform*!
_class
loc:@dense_3/kernel*
dtype0

"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel*!
_class
loc:@dense_3/kernel*
dtype0* 
_output_shapes
:

\
dense_3/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
Ѓ
dense_3/biasVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_3/bias*
_class
loc:@dense_3/bias*
	container *
shape:
i
-dense_3/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_3/bias*
_output_shapes
: 
r
dense_3/bias/AssignAssignVariableOpdense_3/biasdense_3/Const*
_class
loc:@dense_3/bias*
dtype0

 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
dtype0*
_output_shapes	
:*
_class
loc:@dense_3/bias
n
dense_3/MatMul/ReadVariableOpReadVariableOpdense_3/kernel*
dtype0* 
_output_shapes
:

Ѓ
dense_3/MatMulMatMulactivation_2/Reludense_3/MatMul/ReadVariableOp*(
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( *
T0
h
dense_3/BiasAdd/ReadVariableOpReadVariableOpdense_3/bias*
dtype0*
_output_shapes	
:

dense_3/BiasAddBiasAdddense_3/MatMuldense_3/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ
n
dropout_3/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_3/cond/switch_tIdentitydropout_3/cond/Switch:1*
T0
*
_output_shapes
: 
[
dropout_3/cond/switch_fIdentitydropout_3/cond/Switch*
T0
*
_output_shapes
: 
Y
dropout_3/cond/pred_idIdentitykeras_learning_phase*
T0
*
_output_shapes
: 
z
dropout_3/cond/dropout/rateConst^dropout_3/cond/switch_t*
valueB
 *ЭЬЬ=*
dtype0*
_output_shapes
: 

dropout_3/cond/dropout/ShapeShape%dropout_3/cond/dropout/Shape/Switch:1*
T0*
out_type0*
_output_shapes
:
С
#dropout_3/cond/dropout/Shape/SwitchSwitchdense_3/BiasAdddropout_3/cond/pred_id*
T0*"
_class
loc:@dense_3/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

)dropout_3/cond/dropout/random_uniform/minConst^dropout_3/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_3/cond/dropout/random_uniform/maxConst^dropout_3/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *  ?
С
3dropout_3/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_3/cond/dropout/Shape*
T0*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2ш­у*
seedБџх)
Ї
)dropout_3/cond/dropout/random_uniform/subSub)dropout_3/cond/dropout/random_uniform/max)dropout_3/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
У
)dropout_3/cond/dropout/random_uniform/mulMul3dropout_3/cond/dropout/random_uniform/RandomUniform)dropout_3/cond/dropout/random_uniform/sub*
T0*(
_output_shapes
:џџџџџџџџџ
Е
%dropout_3/cond/dropout/random_uniformAdd)dropout_3/cond/dropout/random_uniform/mul)dropout_3/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_3/cond/dropout/sub/xConst^dropout_3/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_3/cond/dropout/subSubdropout_3/cond/dropout/sub/xdropout_3/cond/dropout/rate*
_output_shapes
: *
T0

 dropout_3/cond/dropout/truediv/xConst^dropout_3/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

dropout_3/cond/dropout/truedivRealDiv dropout_3/cond/dropout/truediv/xdropout_3/cond/dropout/sub*
T0*
_output_shapes
: 
Њ
#dropout_3/cond/dropout/GreaterEqualGreaterEqual%dropout_3/cond/dropout/random_uniformdropout_3/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_3/cond/dropout/mulMul%dropout_3/cond/dropout/Shape/Switch:1dropout_3/cond/dropout/truediv*(
_output_shapes
:џџџџџџџџџ*
T0

dropout_3/cond/dropout/CastCast#dropout_3/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0

dropout_3/cond/dropout/mul_1Muldropout_3/cond/dropout/muldropout_3/cond/dropout/Cast*
T0*(
_output_shapes
:џџџџџџџџџ
Е
dropout_3/cond/Switch_1Switchdense_3/BiasAdddropout_3/cond/pred_id*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ*
T0*"
_class
loc:@dense_3/BiasAdd

dropout_3/cond/MergeMergedropout_3/cond/Switch_1dropout_3/cond/dropout/mul_1*
T0*
N**
_output_shapes
:џџџџџџџџџ: 
b
activation_3/ReluReludropout_3/cond/Merge*
T0*(
_output_shapes
:џџџџџџџџџ
m
dense_4/random_uniform/shapeConst*
dtype0*
_output_shapes
:*
valueB"      
_
dense_4/random_uniform/minConst*
valueB
 *   О*
dtype0*
_output_shapes
: 
_
dense_4/random_uniform/maxConst*
valueB
 *   >*
dtype0*
_output_shapes
: 
Љ
$dense_4/random_uniform/RandomUniformRandomUniformdense_4/random_uniform/shape*
seedБџх)*
T0*
dtype0* 
_output_shapes
:
*
seed2Ъ
z
dense_4/random_uniform/subSubdense_4/random_uniform/maxdense_4/random_uniform/min*
T0*
_output_shapes
: 

dense_4/random_uniform/mulMul$dense_4/random_uniform/RandomUniformdense_4/random_uniform/sub*
T0* 
_output_shapes
:


dense_4/random_uniformAdddense_4/random_uniform/muldense_4/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_4/kernelVarHandleOp*
shared_namedense_4/kernel*!
_class
loc:@dense_4/kernel*
	container *
shape:
*
dtype0*
_output_shapes
: 
m
/dense_4/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_4/kernel*
_output_shapes
: 

dense_4/kernel/AssignAssignVariableOpdense_4/kerneldense_4/random_uniform*!
_class
loc:@dense_4/kernel*
dtype0

"dense_4/kernel/Read/ReadVariableOpReadVariableOpdense_4/kernel*!
_class
loc:@dense_4/kernel*
dtype0* 
_output_shapes
:

\
dense_4/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
Ѓ
dense_4/biasVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_4/bias*
_class
loc:@dense_4/bias*
	container *
shape:
i
-dense_4/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_4/bias*
_output_shapes
: 
r
dense_4/bias/AssignAssignVariableOpdense_4/biasdense_4/Const*
dtype0*
_class
loc:@dense_4/bias

 dense_4/bias/Read/ReadVariableOpReadVariableOpdense_4/bias*
_class
loc:@dense_4/bias*
dtype0*
_output_shapes	
:
n
dense_4/MatMul/ReadVariableOpReadVariableOpdense_4/kernel*
dtype0* 
_output_shapes
:

Ѓ
dense_4/MatMulMatMulactivation_3/Reludense_4/MatMul/ReadVariableOp*
transpose_b( *
T0*(
_output_shapes
:џџџџџџџџџ*
transpose_a( 
h
dense_4/BiasAdd/ReadVariableOpReadVariableOpdense_4/bias*
dtype0*
_output_shapes	
:

dense_4/BiasAddBiasAdddense_4/MatMuldense_4/BiasAdd/ReadVariableOp*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ*
T0
]
activation_4/TanhTanhdense_4/BiasAdd*
T0*(
_output_shapes
:џџџџџџџџџ
l
lambda_1/l2_normalize/SquareSquareactivation_4/Tanh*(
_output_shapes
:џџџџџџџџџ*
T0
v
+lambda_1/l2_normalize/Sum/reduction_indicesConst*
valueB :
џџџџџџџџџ*
dtype0*
_output_shapes
: 
К
lambda_1/l2_normalize/SumSumlambda_1/l2_normalize/Square+lambda_1/l2_normalize/Sum/reduction_indices*
	keep_dims(*

Tidx0*
T0*'
_output_shapes
:џџџџџџџџџ
d
lambda_1/l2_normalize/Maximum/yConst*
dtype0*
_output_shapes
: *
valueB
 *ЬМ+

lambda_1/l2_normalize/MaximumMaximumlambda_1/l2_normalize/Sumlambda_1/l2_normalize/Maximum/y*
T0*'
_output_shapes
:џџџџџџџџџ
u
lambda_1/l2_normalize/RsqrtRsqrtlambda_1/l2_normalize/Maximum*
T0*'
_output_shapes
:џџџџџџџџџ

lambda_1/l2_normalizeMulactivation_4/Tanhlambda_1/l2_normalize/Rsqrt*
T0*(
_output_shapes
:џџџџџџџџџ
y
lambda_1/PlaceholderPlaceholder*
dtype0*(
_output_shapes
:џџџџџџџџџ*
shape:џџџџџџџџџ
q
lambda_1/l2_normalize_1/SquareSquarelambda_1/Placeholder*
T0*(
_output_shapes
:џџџџџџџџџ
x
-lambda_1/l2_normalize_1/Sum/reduction_indicesConst*
valueB :
џџџџџџџџџ*
dtype0*
_output_shapes
: 
Р
lambda_1/l2_normalize_1/SumSumlambda_1/l2_normalize_1/Square-lambda_1/l2_normalize_1/Sum/reduction_indices*
T0*'
_output_shapes
:џџџџџџџџџ*
	keep_dims(*

Tidx0
f
!lambda_1/l2_normalize_1/Maximum/yConst*
valueB
 *ЬМ+*
dtype0*
_output_shapes
: 

lambda_1/l2_normalize_1/MaximumMaximumlambda_1/l2_normalize_1/Sum!lambda_1/l2_normalize_1/Maximum/y*
T0*'
_output_shapes
:џџџџџџџџџ
y
lambda_1/l2_normalize_1/RsqrtRsqrtlambda_1/l2_normalize_1/Maximum*
T0*'
_output_shapes
:џџџџџџџџџ

lambda_1/l2_normalize_1Mullambda_1/Placeholderlambda_1/l2_normalize_1/Rsqrt*
T0*(
_output_shapes
:џџџџџџџџџ

)Adam/iterations/Initializer/initial_valueConst*
value	B	 R *"
_class
loc:@Adam/iterations*
dtype0	*
_output_shapes
: 
Ї
Adam/iterationsVarHandleOp*"
_class
loc:@Adam/iterations*
	container *
shape: *
dtype0	*
_output_shapes
: * 
shared_nameAdam/iterations
o
0Adam/iterations/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam/iterations*
_output_shapes
: 

Adam/iterations/AssignAssignVariableOpAdam/iterations)Adam/iterations/Initializer/initial_value*"
_class
loc:@Adam/iterations*
dtype0	

#Adam/iterations/Read/ReadVariableOpReadVariableOpAdam/iterations*
dtype0	*
_output_shapes
: *"
_class
loc:@Adam/iterations

,Adam/learning_rate/Initializer/initial_valueConst*
dtype0*
_output_shapes
: *
valueB
 *o:*%
_class
loc:@Adam/learning_rate
А
Adam/learning_rateVarHandleOp*
dtype0*
_output_shapes
: *#
shared_nameAdam/learning_rate*%
_class
loc:@Adam/learning_rate*
	container *
shape: 
u
3Adam/learning_rate/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam/learning_rate*
_output_shapes
: 
Ѓ
Adam/learning_rate/AssignAssignVariableOpAdam/learning_rate,Adam/learning_rate/Initializer/initial_value*
dtype0*%
_class
loc:@Adam/learning_rate

&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*%
_class
loc:@Adam/learning_rate*
dtype0*
_output_shapes
: 

%Adam/beta_1/Initializer/initial_valueConst*
dtype0*
_output_shapes
: *
valueB
 *fff?*
_class
loc:@Adam/beta_1

Adam/beta_1VarHandleOp*
dtype0*
_output_shapes
: *
shared_nameAdam/beta_1*
_class
loc:@Adam/beta_1*
	container *
shape: 
g
,Adam/beta_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam/beta_1*
_output_shapes
: 

Adam/beta_1/AssignAssignVariableOpAdam/beta_1%Adam/beta_1/Initializer/initial_value*
_class
loc:@Adam/beta_1*
dtype0

Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
dtype0*
_output_shapes
: *
_class
loc:@Adam/beta_1

%Adam/beta_2/Initializer/initial_valueConst*
dtype0*
_output_shapes
: *
valueB
 *wО?*
_class
loc:@Adam/beta_2

Adam/beta_2VarHandleOp*
	container *
shape: *
dtype0*
_output_shapes
: *
shared_nameAdam/beta_2*
_class
loc:@Adam/beta_2
g
,Adam/beta_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam/beta_2*
_output_shapes
: 

Adam/beta_2/AssignAssignVariableOpAdam/beta_2%Adam/beta_2/Initializer/initial_value*
dtype0*
_class
loc:@Adam/beta_2

Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_class
loc:@Adam/beta_2*
dtype0*
_output_shapes
: 

$Adam/decay/Initializer/initial_valueConst*
valueB
 *    *
_class
loc:@Adam/decay*
dtype0*
_output_shapes
: 


Adam/decayVarHandleOp*
dtype0*
_output_shapes
: *
shared_name
Adam/decay*
_class
loc:@Adam/decay*
	container *
shape: 
e
+Adam/decay/IsInitialized/VarIsInitializedOpVarIsInitializedOp
Adam/decay*
_output_shapes
: 

Adam/decay/AssignAssignVariableOp
Adam/decay$Adam/decay/Initializer/initial_value*
_class
loc:@Adam/decay*
dtype0

Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_class
loc:@Adam/decay*
dtype0*
_output_shapes
: 

lambda_1_targetPlaceholder*%
shape:џџџџџџџџџџџџџџџџџџ*
dtype0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
r
lambda_1_sample_weightsPlaceholder*
dtype0*#
_output_shapes
:џџџџџџџџџ*
shape:џџџџџџџџџ
J
ConstConst*
valueB
 *    *
dtype0*
_output_shapes
: 

totalVarHandleOp*
	container *
shape: *
dtype0*
_output_shapes
: *
shared_nametotal*
_class

loc:@total
[
&total/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal*
_output_shapes
: 
U
total/AssignAssignVariableOptotalConst*
_class

loc:@total*
dtype0
q
total/Read/ReadVariableOpReadVariableOptotal*
_class

loc:@total*
dtype0*
_output_shapes
: 
L
Const_1Const*
dtype0*
_output_shapes
: *
valueB
 *    

countVarHandleOp*
shared_namecount*
_class

loc:@count*
	container *
shape: *
dtype0*
_output_shapes
: 
[
&count/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount*
_output_shapes
: 
W
count/AssignAssignVariableOpcountConst_1*
dtype0*
_class

loc:@count
q
count/Read/ReadVariableOpReadVariableOpcount*
_class

loc:@count*
dtype0*
_output_shapes
: 
e
#metrics/corr/Mean/reduction_indicesConst*
value	B : *
dtype0*
_output_shapes
: 

metrics/corr/MeanMeanlambda_1_target#metrics/corr/Mean/reduction_indices*
	keep_dims( *

Tidx0*
T0*#
_output_shapes
:џџџџџџџџџ
g
%metrics/corr/Mean_1/reduction_indicesConst*
dtype0*
_output_shapes
: *
value	B : 

metrics/corr/Mean_1Meanlambda_1/l2_normalize%metrics/corr/Mean_1/reduction_indices*
	keep_dims( *

Tidx0*
T0*
_output_shapes	
:
v
metrics/corr/subSublambda_1_targetmetrics/corr/Mean*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
x
metrics/corr/sub_1Sublambda_1/l2_normalizemetrics/corr/Mean_1*(
_output_shapes
:џџџџџџџџџ*
T0
p
metrics/corr/mulMulmetrics/corr/submetrics/corr/sub_1*
T0*(
_output_shapes
:џџџџџџџџџ
c
metrics/corr/ConstConst*
valueB"       *
dtype0*
_output_shapes
:
{
metrics/corr/SumSummetrics/corr/mulmetrics/corr/Const*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
x
metrics/corr/mul_1Mulmetrics/corr/submetrics/corr/sub*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ*
T0
e
metrics/corr/Const_1Const*
valueB"       *
dtype0*
_output_shapes
:

metrics/corr/Sum_1Summetrics/corr/mul_1metrics/corr/Const_1*
_output_shapes
: *
	keep_dims( *

Tidx0*
T0
t
metrics/corr/mul_2Mulmetrics/corr/sub_1metrics/corr/sub_1*
T0*(
_output_shapes
:џџџџџџџџџ
e
metrics/corr/Const_2Const*
valueB"       *
dtype0*
_output_shapes
:

metrics/corr/Sum_2Summetrics/corr/mul_2metrics/corr/Const_2*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
b
metrics/corr/mul_3Mulmetrics/corr/Sum_1metrics/corr/Sum_2*
_output_shapes
: *
T0
Y
metrics/corr/Const_3Const*
dtype0*
_output_shapes
: *
valueB
 *    
Y
metrics/corr/Const_4Const*
valueB
 *  *
dtype0*
_output_shapes
: 
x
"metrics/corr/clip_by_value/MinimumMinimummetrics/corr/mul_3metrics/corr/Const_4*
T0*
_output_shapes
: 

metrics/corr/clip_by_valueMaximum"metrics/corr/clip_by_value/Minimummetrics/corr/Const_3*
T0*
_output_shapes
: 
V
metrics/corr/SqrtSqrtmetrics/corr/clip_by_value*
_output_shapes
: *
T0
e
metrics/corr/truedivRealDivmetrics/corr/Summetrics/corr/Sqrt*
T0*
_output_shapes
: 
W
metrics/corr/Const_5Const*
valueB *
dtype0*
_output_shapes
: 

metrics/corr/Mean_2Meanmetrics/corr/truedivmetrics/corr/Const_5*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
W
metrics/corr/Const_6Const*
valueB *
dtype0*
_output_shapes
: 

metrics/corr/Sum_3Summetrics/corr/Mean_2metrics/corr/Const_6*
	keep_dims( *

Tidx0*
T0*
_output_shapes
: 
_
 metrics/corr/AssignAddVariableOpAssignAddVariableOptotalmetrics/corr/Sum_3*
dtype0
|
metrics/corr/ReadVariableOpReadVariableOptotal!^metrics/corr/AssignAddVariableOp*
dtype0*
_output_shapes
: 
S
metrics/corr/SizeConst*
value	B :*
dtype0*
_output_shapes
: 
l
metrics/corr/CastCastmetrics/corr/Size*
Truncate( *
_output_shapes
: *

DstT0*

SrcT0
`
"metrics/corr/AssignAddVariableOp_1AssignAddVariableOpcountmetrics/corr/Cast*
dtype0

metrics/corr/ReadVariableOp_1ReadVariableOpcount#^metrics/corr/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ѓ
metrics/corr/ReadVariableOp_2ReadVariableOptotal!^metrics/corr/AssignAddVariableOp#^metrics/corr/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ћ
%metrics/corr/truediv_1/ReadVariableOpReadVariableOpcount!^metrics/corr/AssignAddVariableOp#^metrics/corr/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

metrics/corr/truediv_1RealDivmetrics/corr/ReadVariableOp_2%metrics/corr/truediv_1/ReadVariableOp*
T0*
_output_shapes
: 
Z
metrics/corr/IdentityIdentitymetrics/corr/truediv_1*
T0*
_output_shapes
: 

)loss/lambda_1_loss/mean_squared_error/subSublambda_1/l2_normalizelambda_1_target*(
_output_shapes
:џџџџџџџџџ*
T0

,loss/lambda_1_loss/mean_squared_error/SquareSquare)loss/lambda_1_loss/mean_squared_error/sub*
T0*(
_output_shapes
:џџџџџџџџџ

<loss/lambda_1_loss/mean_squared_error/Mean/reduction_indicesConst*
dtype0*
_output_shapes
: *
valueB :
џџџџџџџџџ
щ
*loss/lambda_1_loss/mean_squared_error/MeanMean,loss/lambda_1_loss/mean_squared_error/Square<loss/lambda_1_loss/mean_squared_error/Mean/reduction_indices*#
_output_shapes
:џџџџџџџџџ*
	keep_dims( *

Tidx0*
T0
Б
7loss/lambda_1_loss/mean_squared_error/weighted_loss/mulMullambda_1_sample_weights*loss/lambda_1_loss/mean_squared_error/Mean*
T0*#
_output_shapes
:џџџџџџџџџ

9loss/lambda_1_loss/mean_squared_error/weighted_loss/ConstConst*
valueB: *
dtype0*
_output_shapes
:
№
7loss/lambda_1_loss/mean_squared_error/weighted_loss/SumSum7loss/lambda_1_loss/mean_squared_error/weighted_loss/mul9loss/lambda_1_loss/mean_squared_error/weighted_loss/Const*
	keep_dims( *

Tidx0*
T0*
_output_shapes
: 
З
Eloss/lambda_1_loss/mean_squared_error/weighted_loss/num_elements/SizeSize7loss/lambda_1_loss/mean_squared_error/weighted_loss/mul*
T0*
out_type0*
_output_shapes
: 
д
Eloss/lambda_1_loss/mean_squared_error/weighted_loss/num_elements/CastCastEloss/lambda_1_loss/mean_squared_error/weighted_loss/num_elements/Size*

SrcT0*
Truncate( *
_output_shapes
: *

DstT0
ч
;loss/lambda_1_loss/mean_squared_error/weighted_loss/truedivRealDiv7loss/lambda_1_loss/mean_squared_error/weighted_loss/SumEloss/lambda_1_loss/mean_squared_error/weighted_loss/num_elements/Cast*
T0*
_output_shapes
: 
O

loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
y
loss/mulMul
loss/mul/x;loss/lambda_1_loss/mean_squared_error/weighted_loss/truediv*
T0*
_output_shapes
: 
J
Const_2Const*
valueB *
dtype0*
_output_shapes
: 
]
MeanMeanloss/mulConst_2*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
`
PlaceholderPlaceholder*
dtype0* 
_output_shapes
:
*
shape:

N
AssignVariableOpAssignVariableOpdense_1/kernelPlaceholder*
dtype0
r
ReadVariableOpReadVariableOpdense_1/kernel^AssignVariableOp*
dtype0* 
_output_shapes
:

X
Placeholder_1Placeholder*
dtype0*
_output_shapes	
:*
shape:
P
AssignVariableOp_1AssignVariableOpdense_1/biasPlaceholder_1*
dtype0
o
ReadVariableOp_1ReadVariableOpdense_1/bias^AssignVariableOp_1*
dtype0*
_output_shapes	
:
b
Placeholder_2Placeholder*
dtype0* 
_output_shapes
:
*
shape:

R
AssignVariableOp_2AssignVariableOpdense_2/kernelPlaceholder_2*
dtype0
v
ReadVariableOp_2ReadVariableOpdense_2/kernel^AssignVariableOp_2*
dtype0* 
_output_shapes
:

X
Placeholder_3Placeholder*
shape:*
dtype0*
_output_shapes	
:
P
AssignVariableOp_3AssignVariableOpdense_2/biasPlaceholder_3*
dtype0
o
ReadVariableOp_3ReadVariableOpdense_2/bias^AssignVariableOp_3*
dtype0*
_output_shapes	
:
b
Placeholder_4Placeholder*
shape:
*
dtype0* 
_output_shapes
:

R
AssignVariableOp_4AssignVariableOpdense_3/kernelPlaceholder_4*
dtype0
v
ReadVariableOp_4ReadVariableOpdense_3/kernel^AssignVariableOp_4*
dtype0* 
_output_shapes
:

X
Placeholder_5Placeholder*
dtype0*
_output_shapes	
:*
shape:
P
AssignVariableOp_5AssignVariableOpdense_3/biasPlaceholder_5*
dtype0
o
ReadVariableOp_5ReadVariableOpdense_3/bias^AssignVariableOp_5*
dtype0*
_output_shapes	
:
b
Placeholder_6Placeholder*
dtype0* 
_output_shapes
:
*
shape:

R
AssignVariableOp_6AssignVariableOpdense_4/kernelPlaceholder_6*
dtype0
v
ReadVariableOp_6ReadVariableOpdense_4/kernel^AssignVariableOp_6*
dtype0* 
_output_shapes
:

X
Placeholder_7Placeholder*
dtype0*
_output_shapes	
:*
shape:
P
AssignVariableOp_7AssignVariableOpdense_4/biasPlaceholder_7*
dtype0
o
ReadVariableOp_7ReadVariableOpdense_4/bias^AssignVariableOp_7*
dtype0*
_output_shapes	
:
N
VarIsInitializedOpVarIsInitializedOpdense_1/bias*
_output_shapes
: 
P
VarIsInitializedOp_1VarIsInitializedOpdense_3/bias*
_output_shapes
: 
O
VarIsInitializedOp_2VarIsInitializedOpAdam/beta_1*
_output_shapes
: 
R
VarIsInitializedOp_3VarIsInitializedOpdense_2/kernel*
_output_shapes
: 
V
VarIsInitializedOp_4VarIsInitializedOpAdam/learning_rate*
_output_shapes
: 
I
VarIsInitializedOp_5VarIsInitializedOptotal*
_output_shapes
: 
R
VarIsInitializedOp_6VarIsInitializedOpdense_1/kernel*
_output_shapes
: 
P
VarIsInitializedOp_7VarIsInitializedOpdense_4/bias*
_output_shapes
: 
R
VarIsInitializedOp_8VarIsInitializedOpdense_4/kernel*
_output_shapes
: 
O
VarIsInitializedOp_9VarIsInitializedOpAdam/beta_2*
_output_shapes
: 
S
VarIsInitializedOp_10VarIsInitializedOpdense_3/kernel*
_output_shapes
: 
Q
VarIsInitializedOp_11VarIsInitializedOpdense_2/bias*
_output_shapes
: 
O
VarIsInitializedOp_12VarIsInitializedOp
Adam/decay*
_output_shapes
: 
T
VarIsInitializedOp_13VarIsInitializedOpAdam/iterations*
_output_shapes
: 
J
VarIsInitializedOp_14VarIsInitializedOpcount*
_output_shapes
: 
е
initNoOp^Adam/beta_1/Assign^Adam/beta_2/Assign^Adam/decay/Assign^Adam/iterations/Assign^Adam/learning_rate/Assign^count/Assign^dense_1/bias/Assign^dense_1/kernel/Assign^dense_2/bias/Assign^dense_2/kernel/Assign^dense_3/bias/Assign^dense_3/kernel/Assign^dense_4/bias/Assign^dense_4/kernel/Assign^total/Assign
t
dropout_4_inputPlaceholder*
dtype0*(
_output_shapes
:џџџџџџџџџ*
shape:џџџџџџџџџ
n
dropout_4/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_4/cond/switch_tIdentitydropout_4/cond/Switch:1*
T0
*
_output_shapes
: 
[
dropout_4/cond/switch_fIdentitydropout_4/cond/Switch*
T0
*
_output_shapes
: 
Y
dropout_4/cond/pred_idIdentitykeras_learning_phase*
T0
*
_output_shapes
: 
z
dropout_4/cond/dropout/rateConst^dropout_4/cond/switch_t*
valueB
 *   ?*
dtype0*
_output_shapes
: 

dropout_4/cond/dropout/ShapeShape%dropout_4/cond/dropout/Shape/Switch:1*
_output_shapes
:*
T0*
out_type0
С
#dropout_4/cond/dropout/Shape/SwitchSwitchdropout_4_inputdropout_4/cond/pred_id*
T0*"
_class
loc:@dropout_4_input*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

)dropout_4/cond/dropout/random_uniform/minConst^dropout_4/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_4/cond/dropout/random_uniform/maxConst^dropout_4/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
С
3dropout_4/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_4/cond/dropout/Shape*
T0*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2ХХ*
seedБџх)
Ї
)dropout_4/cond/dropout/random_uniform/subSub)dropout_4/cond/dropout/random_uniform/max)dropout_4/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
У
)dropout_4/cond/dropout/random_uniform/mulMul3dropout_4/cond/dropout/random_uniform/RandomUniform)dropout_4/cond/dropout/random_uniform/sub*
T0*(
_output_shapes
:џџџџџџџџџ
Е
%dropout_4/cond/dropout/random_uniformAdd)dropout_4/cond/dropout/random_uniform/mul)dropout_4/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_4/cond/dropout/sub/xConst^dropout_4/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_4/cond/dropout/subSubdropout_4/cond/dropout/sub/xdropout_4/cond/dropout/rate*
T0*
_output_shapes
: 

 dropout_4/cond/dropout/truediv/xConst^dropout_4/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *  ?

dropout_4/cond/dropout/truedivRealDiv dropout_4/cond/dropout/truediv/xdropout_4/cond/dropout/sub*
_output_shapes
: *
T0
Њ
#dropout_4/cond/dropout/GreaterEqualGreaterEqual%dropout_4/cond/dropout/random_uniformdropout_4/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_4/cond/dropout/mulMul%dropout_4/cond/dropout/Shape/Switch:1dropout_4/cond/dropout/truediv*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_4/cond/dropout/CastCast#dropout_4/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0

dropout_4/cond/dropout/mul_1Muldropout_4/cond/dropout/muldropout_4/cond/dropout/Cast*
T0*(
_output_shapes
:џџџџџџџџџ
Е
dropout_4/cond/Switch_1Switchdropout_4_inputdropout_4/cond/pred_id*
T0*"
_class
loc:@dropout_4_input*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

dropout_4/cond/MergeMergedropout_4/cond/Switch_1dropout_4/cond/dropout/mul_1*
T0*
N**
_output_shapes
:џџџџџџџџџ: 
m
dense_5/random_uniform/shapeConst*
dtype0*
_output_shapes
:*
valueB"      
_
dense_5/random_uniform/minConst*
valueB
 *ѓ5Н*
dtype0*
_output_shapes
: 
_
dense_5/random_uniform/maxConst*
valueB
 *ѓ5=*
dtype0*
_output_shapes
: 
Љ
$dense_5/random_uniform/RandomUniformRandomUniformdense_5/random_uniform/shape*
T0*
dtype0* 
_output_shapes
:
*
seed2НЂw*
seedБџх)
z
dense_5/random_uniform/subSubdense_5/random_uniform/maxdense_5/random_uniform/min*
T0*
_output_shapes
: 

dense_5/random_uniform/mulMul$dense_5/random_uniform/RandomUniformdense_5/random_uniform/sub*
T0* 
_output_shapes
:


dense_5/random_uniformAdddense_5/random_uniform/muldense_5/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_5/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_5/kernel*!
_class
loc:@dense_5/kernel*
	container *
shape:

m
/dense_5/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_5/kernel*
_output_shapes
: 

dense_5/kernel/AssignAssignVariableOpdense_5/kerneldense_5/random_uniform*!
_class
loc:@dense_5/kernel*
dtype0

"dense_5/kernel/Read/ReadVariableOpReadVariableOpdense_5/kernel*!
_class
loc:@dense_5/kernel*
dtype0* 
_output_shapes
:

\
dense_5/ConstConst*
dtype0*
_output_shapes	
:*
valueB*    
Ѓ
dense_5/biasVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_5/bias*
_class
loc:@dense_5/bias*
	container *
shape:
i
-dense_5/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_5/bias*
_output_shapes
: 
r
dense_5/bias/AssignAssignVariableOpdense_5/biasdense_5/Const*
_class
loc:@dense_5/bias*
dtype0

 dense_5/bias/Read/ReadVariableOpReadVariableOpdense_5/bias*
dtype0*
_output_shapes	
:*
_class
loc:@dense_5/bias
n
dense_5/MatMul/ReadVariableOpReadVariableOpdense_5/kernel*
dtype0* 
_output_shapes
:

І
dense_5/MatMulMatMuldropout_4/cond/Mergedense_5/MatMul/ReadVariableOp*
T0*(
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( 
h
dense_5/BiasAdd/ReadVariableOpReadVariableOpdense_5/bias*
dtype0*
_output_shapes	
:

dense_5/BiasAddBiasAdddense_5/MatMuldense_5/BiasAdd/ReadVariableOp*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ*
T0
n
dropout_5/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_5/cond/switch_tIdentitydropout_5/cond/Switch:1*
_output_shapes
: *
T0

[
dropout_5/cond/switch_fIdentitydropout_5/cond/Switch*
T0
*
_output_shapes
: 
Y
dropout_5/cond/pred_idIdentitykeras_learning_phase*
T0
*
_output_shapes
: 
z
dropout_5/cond/dropout/rateConst^dropout_5/cond/switch_t*
valueB
 *ЭЬL>*
dtype0*
_output_shapes
: 

dropout_5/cond/dropout/ShapeShape%dropout_5/cond/dropout/Shape/Switch:1*
T0*
out_type0*
_output_shapes
:
С
#dropout_5/cond/dropout/Shape/SwitchSwitchdense_5/BiasAdddropout_5/cond/pred_id*
T0*"
_class
loc:@dense_5/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

)dropout_5/cond/dropout/random_uniform/minConst^dropout_5/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_5/cond/dropout/random_uniform/maxConst^dropout_5/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *  ?
С
3dropout_5/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_5/cond/dropout/Shape*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2рљ*
seedБџх)*
T0
Ї
)dropout_5/cond/dropout/random_uniform/subSub)dropout_5/cond/dropout/random_uniform/max)dropout_5/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
У
)dropout_5/cond/dropout/random_uniform/mulMul3dropout_5/cond/dropout/random_uniform/RandomUniform)dropout_5/cond/dropout/random_uniform/sub*(
_output_shapes
:џџџџџџџџџ*
T0
Е
%dropout_5/cond/dropout/random_uniformAdd)dropout_5/cond/dropout/random_uniform/mul)dropout_5/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_5/cond/dropout/sub/xConst^dropout_5/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_5/cond/dropout/subSubdropout_5/cond/dropout/sub/xdropout_5/cond/dropout/rate*
_output_shapes
: *
T0

 dropout_5/cond/dropout/truediv/xConst^dropout_5/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *  ?

dropout_5/cond/dropout/truedivRealDiv dropout_5/cond/dropout/truediv/xdropout_5/cond/dropout/sub*
T0*
_output_shapes
: 
Њ
#dropout_5/cond/dropout/GreaterEqualGreaterEqual%dropout_5/cond/dropout/random_uniformdropout_5/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_5/cond/dropout/mulMul%dropout_5/cond/dropout/Shape/Switch:1dropout_5/cond/dropout/truediv*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_5/cond/dropout/CastCast#dropout_5/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0

dropout_5/cond/dropout/mul_1Muldropout_5/cond/dropout/muldropout_5/cond/dropout/Cast*(
_output_shapes
:џџџџџџџџџ*
T0
Е
dropout_5/cond/Switch_1Switchdense_5/BiasAdddropout_5/cond/pred_id*
T0*"
_class
loc:@dense_5/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

dropout_5/cond/MergeMergedropout_5/cond/Switch_1dropout_5/cond/dropout/mul_1*
T0*
N**
_output_shapes
:џџџџџџџџџ: 
b
activation_5/ReluReludropout_5/cond/Merge*(
_output_shapes
:џџџџџџџџџ*
T0
m
dense_6/random_uniform/shapeConst*
valueB"      *
dtype0*
_output_shapes
:
_
dense_6/random_uniform/minConst*
valueB
 *  Н*
dtype0*
_output_shapes
: 
_
dense_6/random_uniform/maxConst*
valueB
 *  =*
dtype0*
_output_shapes
: 
Њ
$dense_6/random_uniform/RandomUniformRandomUniformdense_6/random_uniform/shape*
T0*
dtype0* 
_output_shapes
:
*
seed2ЋхЗ*
seedБџх)
z
dense_6/random_uniform/subSubdense_6/random_uniform/maxdense_6/random_uniform/min*
_output_shapes
: *
T0

dense_6/random_uniform/mulMul$dense_6/random_uniform/RandomUniformdense_6/random_uniform/sub*
T0* 
_output_shapes
:


dense_6/random_uniformAdddense_6/random_uniform/muldense_6/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_6/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_6/kernel*!
_class
loc:@dense_6/kernel*
	container *
shape:

m
/dense_6/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_6/kernel*
_output_shapes
: 

dense_6/kernel/AssignAssignVariableOpdense_6/kerneldense_6/random_uniform*!
_class
loc:@dense_6/kernel*
dtype0

"dense_6/kernel/Read/ReadVariableOpReadVariableOpdense_6/kernel*
dtype0* 
_output_shapes
:
*!
_class
loc:@dense_6/kernel
\
dense_6/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
Ѓ
dense_6/biasVarHandleOp*
shared_namedense_6/bias*
_class
loc:@dense_6/bias*
	container *
shape:*
dtype0*
_output_shapes
: 
i
-dense_6/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_6/bias*
_output_shapes
: 
r
dense_6/bias/AssignAssignVariableOpdense_6/biasdense_6/Const*
dtype0*
_class
loc:@dense_6/bias

 dense_6/bias/Read/ReadVariableOpReadVariableOpdense_6/bias*
_class
loc:@dense_6/bias*
dtype0*
_output_shapes	
:
n
dense_6/MatMul/ReadVariableOpReadVariableOpdense_6/kernel*
dtype0* 
_output_shapes
:

Ѓ
dense_6/MatMulMatMulactivation_5/Reludense_6/MatMul/ReadVariableOp*(
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( *
T0
h
dense_6/BiasAdd/ReadVariableOpReadVariableOpdense_6/bias*
dtype0*
_output_shapes	
:

dense_6/BiasAddBiasAdddense_6/MatMuldense_6/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ
n
dropout_6/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_6/cond/switch_tIdentitydropout_6/cond/Switch:1*
_output_shapes
: *
T0

[
dropout_6/cond/switch_fIdentitydropout_6/cond/Switch*
T0
*
_output_shapes
: 
Y
dropout_6/cond/pred_idIdentitykeras_learning_phase*
T0
*
_output_shapes
: 
z
dropout_6/cond/dropout/rateConst^dropout_6/cond/switch_t*
valueB
 *ЭЬL>*
dtype0*
_output_shapes
: 

dropout_6/cond/dropout/ShapeShape%dropout_6/cond/dropout/Shape/Switch:1*
T0*
out_type0*
_output_shapes
:
С
#dropout_6/cond/dropout/Shape/SwitchSwitchdense_6/BiasAdddropout_6/cond/pred_id*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ*
T0*"
_class
loc:@dense_6/BiasAdd

)dropout_6/cond/dropout/random_uniform/minConst^dropout_6/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_6/cond/dropout/random_uniform/maxConst^dropout_6/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
С
3dropout_6/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_6/cond/dropout/Shape*
seedБџх)*
T0*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2ќ
Ї
)dropout_6/cond/dropout/random_uniform/subSub)dropout_6/cond/dropout/random_uniform/max)dropout_6/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
У
)dropout_6/cond/dropout/random_uniform/mulMul3dropout_6/cond/dropout/random_uniform/RandomUniform)dropout_6/cond/dropout/random_uniform/sub*
T0*(
_output_shapes
:џџџџџџџџџ
Е
%dropout_6/cond/dropout/random_uniformAdd)dropout_6/cond/dropout/random_uniform/mul)dropout_6/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_6/cond/dropout/sub/xConst^dropout_6/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_6/cond/dropout/subSubdropout_6/cond/dropout/sub/xdropout_6/cond/dropout/rate*
T0*
_output_shapes
: 

 dropout_6/cond/dropout/truediv/xConst^dropout_6/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

dropout_6/cond/dropout/truedivRealDiv dropout_6/cond/dropout/truediv/xdropout_6/cond/dropout/sub*
T0*
_output_shapes
: 
Њ
#dropout_6/cond/dropout/GreaterEqualGreaterEqual%dropout_6/cond/dropout/random_uniformdropout_6/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_6/cond/dropout/mulMul%dropout_6/cond/dropout/Shape/Switch:1dropout_6/cond/dropout/truediv*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_6/cond/dropout/CastCast#dropout_6/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0

dropout_6/cond/dropout/mul_1Muldropout_6/cond/dropout/muldropout_6/cond/dropout/Cast*(
_output_shapes
:џџџџџџџџџ*
T0
Е
dropout_6/cond/Switch_1Switchdense_6/BiasAdddropout_6/cond/pred_id*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ*
T0*"
_class
loc:@dense_6/BiasAdd

dropout_6/cond/MergeMergedropout_6/cond/Switch_1dropout_6/cond/dropout/mul_1*
N**
_output_shapes
:џџџџџџџџџ: *
T0
b
activation_6/ReluReludropout_6/cond/Merge*
T0*(
_output_shapes
:џџџџџџџџџ
m
dense_7/random_uniform/shapeConst*
dtype0*
_output_shapes
:*
valueB"      
_
dense_7/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *ѓЕН
_
dense_7/random_uniform/maxConst*
valueB
 *ѓЕ=*
dtype0*
_output_shapes
: 
Њ
$dense_7/random_uniform/RandomUniformRandomUniformdense_7/random_uniform/shape*
T0*
dtype0* 
_output_shapes
:
*
seed2Ю*
seedБџх)
z
dense_7/random_uniform/subSubdense_7/random_uniform/maxdense_7/random_uniform/min*
_output_shapes
: *
T0

dense_7/random_uniform/mulMul$dense_7/random_uniform/RandomUniformdense_7/random_uniform/sub*
T0* 
_output_shapes
:


dense_7/random_uniformAdddense_7/random_uniform/muldense_7/random_uniform/min*
T0* 
_output_shapes
:

Ў
dense_7/kernelVarHandleOp*
shared_namedense_7/kernel*!
_class
loc:@dense_7/kernel*
	container *
shape:
*
dtype0*
_output_shapes
: 
m
/dense_7/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_7/kernel*
_output_shapes
: 

dense_7/kernel/AssignAssignVariableOpdense_7/kerneldense_7/random_uniform*!
_class
loc:@dense_7/kernel*
dtype0

"dense_7/kernel/Read/ReadVariableOpReadVariableOpdense_7/kernel*!
_class
loc:@dense_7/kernel*
dtype0* 
_output_shapes
:

\
dense_7/ConstConst*
dtype0*
_output_shapes	
:*
valueB*    
Ѓ
dense_7/biasVarHandleOp*
_class
loc:@dense_7/bias*
	container *
shape:*
dtype0*
_output_shapes
: *
shared_namedense_7/bias
i
-dense_7/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_7/bias*
_output_shapes
: 
r
dense_7/bias/AssignAssignVariableOpdense_7/biasdense_7/Const*
_class
loc:@dense_7/bias*
dtype0

 dense_7/bias/Read/ReadVariableOpReadVariableOpdense_7/bias*
_class
loc:@dense_7/bias*
dtype0*
_output_shapes	
:
n
dense_7/MatMul/ReadVariableOpReadVariableOpdense_7/kernel*
dtype0* 
_output_shapes
:

Ѓ
dense_7/MatMulMatMulactivation_6/Reludense_7/MatMul/ReadVariableOp*
T0*(
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( 
h
dense_7/BiasAdd/ReadVariableOpReadVariableOpdense_7/bias*
dtype0*
_output_shapes	
:

dense_7/BiasAddBiasAdddense_7/MatMuldense_7/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC*(
_output_shapes
:џџџџџџџџџ
n
dropout_7/cond/SwitchSwitchkeras_learning_phasekeras_learning_phase*
T0
*
_output_shapes
: : 
]
dropout_7/cond/switch_tIdentitydropout_7/cond/Switch:1*
_output_shapes
: *
T0

[
dropout_7/cond/switch_fIdentitydropout_7/cond/Switch*
T0
*
_output_shapes
: 
Y
dropout_7/cond/pred_idIdentitykeras_learning_phase*
_output_shapes
: *
T0

z
dropout_7/cond/dropout/rateConst^dropout_7/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *ЭЬL>

dropout_7/cond/dropout/ShapeShape%dropout_7/cond/dropout/Shape/Switch:1*
T0*
out_type0*
_output_shapes
:
С
#dropout_7/cond/dropout/Shape/SwitchSwitchdense_7/BiasAdddropout_7/cond/pred_id*
T0*"
_class
loc:@dense_7/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

)dropout_7/cond/dropout/random_uniform/minConst^dropout_7/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

)dropout_7/cond/dropout/random_uniform/maxConst^dropout_7/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Р
3dropout_7/cond/dropout/random_uniform/RandomUniformRandomUniformdropout_7/cond/dropout/Shape*
dtype0*(
_output_shapes
:џџџџџџџџџ*
seed2КОT*
seedБџх)*
T0
Ї
)dropout_7/cond/dropout/random_uniform/subSub)dropout_7/cond/dropout/random_uniform/max)dropout_7/cond/dropout/random_uniform/min*
_output_shapes
: *
T0
У
)dropout_7/cond/dropout/random_uniform/mulMul3dropout_7/cond/dropout/random_uniform/RandomUniform)dropout_7/cond/dropout/random_uniform/sub*
T0*(
_output_shapes
:џџџџџџџџџ
Е
%dropout_7/cond/dropout/random_uniformAdd)dropout_7/cond/dropout/random_uniform/mul)dropout_7/cond/dropout/random_uniform/min*
T0*(
_output_shapes
:џџџџџџџџџ
{
dropout_7/cond/dropout/sub/xConst^dropout_7/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
}
dropout_7/cond/dropout/subSubdropout_7/cond/dropout/sub/xdropout_7/cond/dropout/rate*
T0*
_output_shapes
: 

 dropout_7/cond/dropout/truediv/xConst^dropout_7/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

dropout_7/cond/dropout/truedivRealDiv dropout_7/cond/dropout/truediv/xdropout_7/cond/dropout/sub*
T0*
_output_shapes
: 
Њ
#dropout_7/cond/dropout/GreaterEqualGreaterEqual%dropout_7/cond/dropout/random_uniformdropout_7/cond/dropout/rate*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_7/cond/dropout/mulMul%dropout_7/cond/dropout/Shape/Switch:1dropout_7/cond/dropout/truediv*
T0*(
_output_shapes
:џџџџџџџџџ

dropout_7/cond/dropout/CastCast#dropout_7/cond/dropout/GreaterEqual*
Truncate( *(
_output_shapes
:џџџџџџџџџ*

DstT0*

SrcT0


dropout_7/cond/dropout/mul_1Muldropout_7/cond/dropout/muldropout_7/cond/dropout/Cast*(
_output_shapes
:џџџџџџџџџ*
T0
Е
dropout_7/cond/Switch_1Switchdense_7/BiasAdddropout_7/cond/pred_id*
T0*"
_class
loc:@dense_7/BiasAdd*<
_output_shapes*
(:џџџџџџџџџ:џџџџџџџџџ

dropout_7/cond/MergeMergedropout_7/cond/Switch_1dropout_7/cond/dropout/mul_1*
T0*
N**
_output_shapes
:џџџџџџџџџ: 
b
activation_7/ReluReludropout_7/cond/Merge*(
_output_shapes
:џџџџџџџџџ*
T0
m
dense_8/random_uniform/shapeConst*
valueB"      *
dtype0*
_output_shapes
:
_
dense_8/random_uniform/minConst*
valueB
 *IvО*
dtype0*
_output_shapes
: 
_
dense_8/random_uniform/maxConst*
dtype0*
_output_shapes
: *
valueB
 *Iv>
Ј
$dense_8/random_uniform/RandomUniformRandomUniformdense_8/random_uniform/shape*
seedБџх)*
T0*
dtype0*
_output_shapes
:	*
seed2н
z
dense_8/random_uniform/subSubdense_8/random_uniform/maxdense_8/random_uniform/min*
T0*
_output_shapes
: 

dense_8/random_uniform/mulMul$dense_8/random_uniform/RandomUniformdense_8/random_uniform/sub*
T0*
_output_shapes
:	

dense_8/random_uniformAdddense_8/random_uniform/muldense_8/random_uniform/min*
T0*
_output_shapes
:	
­
dense_8/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shared_namedense_8/kernel*!
_class
loc:@dense_8/kernel*
	container *
shape:	
m
/dense_8/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_8/kernel*
_output_shapes
: 

dense_8/kernel/AssignAssignVariableOpdense_8/kerneldense_8/random_uniform*!
_class
loc:@dense_8/kernel*
dtype0

"dense_8/kernel/Read/ReadVariableOpReadVariableOpdense_8/kernel*!
_class
loc:@dense_8/kernel*
dtype0*
_output_shapes
:	
Z
dense_8/ConstConst*
valueB*    *
dtype0*
_output_shapes
:
Ђ
dense_8/biasVarHandleOp*
shape:*
dtype0*
_output_shapes
: *
shared_namedense_8/bias*
_class
loc:@dense_8/bias*
	container 
i
-dense_8/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_8/bias*
_output_shapes
: 
r
dense_8/bias/AssignAssignVariableOpdense_8/biasdense_8/Const*
dtype0*
_class
loc:@dense_8/bias

 dense_8/bias/Read/ReadVariableOpReadVariableOpdense_8/bias*
_class
loc:@dense_8/bias*
dtype0*
_output_shapes
:
m
dense_8/MatMul/ReadVariableOpReadVariableOpdense_8/kernel*
dtype0*
_output_shapes
:	
Ђ
dense_8/MatMulMatMulactivation_7/Reludense_8/MatMul/ReadVariableOp*
T0*'
_output_shapes
:џџџџџџџџџ*
transpose_a( *
transpose_b( 
g
dense_8/BiasAdd/ReadVariableOpReadVariableOpdense_8/bias*
dtype0*
_output_shapes
:

dense_8/BiasAddBiasAdddense_8/MatMuldense_8/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC*'
_output_shapes
:џџџџџџџџџ
d
activation_8/IdentityIdentitydense_8/BiasAdd*
T0*'
_output_shapes
:џџџџџџџџџ

+Adam_1/iterations/Initializer/initial_valueConst*
value	B	 R *$
_class
loc:@Adam_1/iterations*
dtype0	*
_output_shapes
: 
­
Adam_1/iterationsVarHandleOp*
shape: *
dtype0	*
_output_shapes
: *"
shared_nameAdam_1/iterations*$
_class
loc:@Adam_1/iterations*
	container 
s
2Adam_1/iterations/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam_1/iterations*
_output_shapes
: 

Adam_1/iterations/AssignAssignVariableOpAdam_1/iterations+Adam_1/iterations/Initializer/initial_value*
dtype0	*$
_class
loc:@Adam_1/iterations

%Adam_1/iterations/Read/ReadVariableOpReadVariableOpAdam_1/iterations*$
_class
loc:@Adam_1/iterations*
dtype0	*
_output_shapes
: 

.Adam_1/learning_rate/Initializer/initial_valueConst*
valueB
 *o:*'
_class
loc:@Adam_1/learning_rate*
dtype0*
_output_shapes
: 
Ж
Adam_1/learning_rateVarHandleOp*'
_class
loc:@Adam_1/learning_rate*
	container *
shape: *
dtype0*
_output_shapes
: *%
shared_nameAdam_1/learning_rate
y
5Adam_1/learning_rate/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam_1/learning_rate*
_output_shapes
: 
Ћ
Adam_1/learning_rate/AssignAssignVariableOpAdam_1/learning_rate.Adam_1/learning_rate/Initializer/initial_value*'
_class
loc:@Adam_1/learning_rate*
dtype0

(Adam_1/learning_rate/Read/ReadVariableOpReadVariableOpAdam_1/learning_rate*
dtype0*
_output_shapes
: *'
_class
loc:@Adam_1/learning_rate

'Adam_1/beta_1/Initializer/initial_valueConst*
valueB
 *fff?* 
_class
loc:@Adam_1/beta_1*
dtype0*
_output_shapes
: 
Ё
Adam_1/beta_1VarHandleOp*
dtype0*
_output_shapes
: *
shared_nameAdam_1/beta_1* 
_class
loc:@Adam_1/beta_1*
	container *
shape: 
k
.Adam_1/beta_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam_1/beta_1*
_output_shapes
: 

Adam_1/beta_1/AssignAssignVariableOpAdam_1/beta_1'Adam_1/beta_1/Initializer/initial_value* 
_class
loc:@Adam_1/beta_1*
dtype0

!Adam_1/beta_1/Read/ReadVariableOpReadVariableOpAdam_1/beta_1* 
_class
loc:@Adam_1/beta_1*
dtype0*
_output_shapes
: 

'Adam_1/beta_2/Initializer/initial_valueConst*
valueB
 *wО?* 
_class
loc:@Adam_1/beta_2*
dtype0*
_output_shapes
: 
Ё
Adam_1/beta_2VarHandleOp*
dtype0*
_output_shapes
: *
shared_nameAdam_1/beta_2* 
_class
loc:@Adam_1/beta_2*
	container *
shape: 
k
.Adam_1/beta_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam_1/beta_2*
_output_shapes
: 

Adam_1/beta_2/AssignAssignVariableOpAdam_1/beta_2'Adam_1/beta_2/Initializer/initial_value* 
_class
loc:@Adam_1/beta_2*
dtype0

!Adam_1/beta_2/Read/ReadVariableOpReadVariableOpAdam_1/beta_2* 
_class
loc:@Adam_1/beta_2*
dtype0*
_output_shapes
: 

&Adam_1/decay/Initializer/initial_valueConst*
valueB
 *    *
_class
loc:@Adam_1/decay*
dtype0*
_output_shapes
: 

Adam_1/decayVarHandleOp*
_class
loc:@Adam_1/decay*
	container *
shape: *
dtype0*
_output_shapes
: *
shared_nameAdam_1/decay
i
-Adam_1/decay/IsInitialized/VarIsInitializedOpVarIsInitializedOpAdam_1/decay*
_output_shapes
: 

Adam_1/decay/AssignAssignVariableOpAdam_1/decay&Adam_1/decay/Initializer/initial_value*
_class
loc:@Adam_1/decay*
dtype0

 Adam_1/decay/Read/ReadVariableOpReadVariableOpAdam_1/decay*
_class
loc:@Adam_1/decay*
dtype0*
_output_shapes
: 
L
Const_3Const*
valueB
 *    *
dtype0*
_output_shapes
: 

total_1VarHandleOp*
dtype0*
_output_shapes
: *
shared_name	total_1*
_class
loc:@total_1*
	container *
shape: 
_
(total_1/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
]
total_1/AssignAssignVariableOptotal_1Const_3*
_class
loc:@total_1*
dtype0
w
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: 
L
Const_4Const*
valueB
 *    *
dtype0*
_output_shapes
: 

count_1VarHandleOp*
shared_name	count_1*
_class
loc:@count_1*
	container *
shape: *
dtype0*
_output_shapes
: 
_
(count_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_1*
_output_shapes
: 
]
count_1/AssignAssignVariableOpcount_1Const_4*
_class
loc:@count_1*
dtype0
w
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_class
loc:@count_1*
dtype0*
_output_shapes
: 
L
Const_5Const*
valueB
 *    *
dtype0*
_output_shapes
: 

total_2VarHandleOp*
dtype0*
_output_shapes
: *
shared_name	total_2*
_class
loc:@total_2*
	container *
shape: 
_
(total_2/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_2*
_output_shapes
: 
]
total_2/AssignAssignVariableOptotal_2Const_5*
_class
loc:@total_2*
dtype0
w
total_2/Read/ReadVariableOpReadVariableOptotal_2*
dtype0*
_output_shapes
: *
_class
loc:@total_2
L
Const_6Const*
valueB
 *    *
dtype0*
_output_shapes
: 

count_2VarHandleOp*
shared_name	count_2*
_class
loc:@count_2*
	container *
shape: *
dtype0*
_output_shapes
: 
_
(count_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_2*
_output_shapes
: 
]
count_2/AssignAssignVariableOpcount_2Const_6*
_class
loc:@count_2*
dtype0
w
count_2/Read/ReadVariableOpReadVariableOpcount_2*
_class
loc:@count_2*
dtype0*
_output_shapes
: 
L
Const_7Const*
valueB
 *    *
dtype0*
_output_shapes
: 

total_3VarHandleOp*
_class
loc:@total_3*
	container *
shape: *
dtype0*
_output_shapes
: *
shared_name	total_3
_
(total_3/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_3*
_output_shapes
: 
]
total_3/AssignAssignVariableOptotal_3Const_7*
_class
loc:@total_3*
dtype0
w
total_3/Read/ReadVariableOpReadVariableOptotal_3*
_class
loc:@total_3*
dtype0*
_output_shapes
: 
L
Const_8Const*
valueB
 *    *
dtype0*
_output_shapes
: 

count_3VarHandleOp*
_class
loc:@count_3*
	container *
shape: *
dtype0*
_output_shapes
: *
shared_name	count_3
_
(count_3/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_3*
_output_shapes
: 
]
count_3/AssignAssignVariableOpcount_3Const_8*
_class
loc:@count_3*
dtype0
w
count_3/Read/ReadVariableOpReadVariableOpcount_3*
_class
loc:@count_3*
dtype0*
_output_shapes
: 
L
Const_9Const*
valueB
 *    *
dtype0*
_output_shapes
: 

total_4VarHandleOp*
dtype0*
_output_shapes
: *
shared_name	total_4*
_class
loc:@total_4*
	container *
shape: 
_
(total_4/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_4*
_output_shapes
: 
]
total_4/AssignAssignVariableOptotal_4Const_9*
dtype0*
_class
loc:@total_4
w
total_4/Read/ReadVariableOpReadVariableOptotal_4*
_class
loc:@total_4*
dtype0*
_output_shapes
: 
M
Const_10Const*
valueB
 *    *
dtype0*
_output_shapes
: 

count_4VarHandleOp*
shape: *
dtype0*
_output_shapes
: *
shared_name	count_4*
_class
loc:@count_4*
	container 
_
(count_4/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_4*
_output_shapes
: 
^
count_4/AssignAssignVariableOpcount_4Const_10*
_class
loc:@count_4*
dtype0
w
count_4/Read/ReadVariableOpReadVariableOpcount_4*
_class
loc:@count_4*
dtype0*
_output_shapes
: 

activation_8_targetPlaceholder*
dtype0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ*%
shape:џџџџџџџџџџџџџџџџџџ
v
activation_8_sample_weightsPlaceholder*
dtype0*#
_output_shapes
:џџџџџџџџџ*
shape:џџџџџџџџџ
M
Const_11Const*
dtype0*
_output_shapes
: *
valueB
 *    

total_5VarHandleOp*
shape: *
dtype0*
_output_shapes
: *
shared_name	total_5*
_class
loc:@total_5*
	container 
_
(total_5/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_5*
_output_shapes
: 
^
total_5/AssignAssignVariableOptotal_5Const_11*
_class
loc:@total_5*
dtype0
w
total_5/Read/ReadVariableOpReadVariableOptotal_5*
_class
loc:@total_5*
dtype0*
_output_shapes
: 
M
Const_12Const*
dtype0*
_output_shapes
: *
valueB
 *    

count_5VarHandleOp*
	container *
shape: *
dtype0*
_output_shapes
: *
shared_name	count_5*
_class
loc:@count_5
_
(count_5/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_5*
_output_shapes
: 
^
count_5/AssignAssignVariableOpcount_5Const_12*
dtype0*
_class
loc:@count_5
w
count_5/Read/ReadVariableOpReadVariableOpcount_5*
dtype0*
_output_shapes
: *
_class
loc:@count_5

metrics_1/rmse/subSubactivation_8/Identityactivation_8_target*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
n
metrics_1/rmse/SquareSquaremetrics_1/rmse/sub*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
e
metrics_1/rmse/ConstConst*
dtype0*
_output_shapes
:*
valueB"       

metrics_1/rmse/SumSummetrics_1/rmse/Squaremetrics_1/rmse/Const*
_output_shapes
: *
	keep_dims( *

Tidx0*
T0
c
"metrics_1/rmse/AssignAddVariableOpAssignAddVariableOptotal_1metrics_1/rmse/Sum*
dtype0

metrics_1/rmse/ReadVariableOpReadVariableOptotal_1#^metrics_1/rmse/AssignAddVariableOp*
dtype0*
_output_shapes
: 
c
metrics_1/rmse/SizeSizemetrics_1/rmse/Square*
T0*
out_type0*
_output_shapes
: 
p
metrics_1/rmse/CastCastmetrics_1/rmse/Size*

SrcT0*
Truncate( *
_output_shapes
: *

DstT0
f
$metrics_1/rmse/AssignAddVariableOp_1AssignAddVariableOpcount_1metrics_1/rmse/Cast*
dtype0

metrics_1/rmse/ReadVariableOp_1ReadVariableOpcount_1%^metrics_1/rmse/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ћ
metrics_1/rmse/ReadVariableOp_2ReadVariableOptotal_1#^metrics_1/rmse/AssignAddVariableOp%^metrics_1/rmse/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Б
%metrics_1/rmse/truediv/ReadVariableOpReadVariableOpcount_1#^metrics_1/rmse/AssignAddVariableOp%^metrics_1/rmse/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

metrics_1/rmse/truedivRealDivmetrics_1/rmse/ReadVariableOp_2%metrics_1/rmse/truediv/ReadVariableOp*
T0*
_output_shapes
: 
Ї
metrics_1/rmse/Const_1Const#^metrics_1/rmse/AssignAddVariableOp%^metrics_1/rmse/AssignAddVariableOp_1*
valueB
 *    *
dtype0*
_output_shapes
: 
Ї
metrics_1/rmse/Const_2Const#^metrics_1/rmse/AssignAddVariableOp%^metrics_1/rmse/AssignAddVariableOp_1*
valueB
 *  *
dtype0*
_output_shapes
: 

$metrics_1/rmse/clip_by_value/MinimumMinimummetrics_1/rmse/truedivmetrics_1/rmse/Const_2*
T0*
_output_shapes
: 

metrics_1/rmse/clip_by_valueMaximum$metrics_1/rmse/clip_by_value/Minimummetrics_1/rmse/Const_1*
T0*
_output_shapes
: 
Z
metrics_1/rmse/SqrtSqrtmetrics_1/rmse/clip_by_value*
T0*
_output_shapes
: 
Y
metrics_1/rmse/IdentityIdentitymetrics_1/rmse/Sqrt*
T0*
_output_shapes
: 

metrics_1/mea/subSubactivation_8/Identityactivation_8_target*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
f
metrics_1/mea/AbsAbsmetrics_1/mea/sub*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
o
$metrics_1/mea/Mean/reduction_indicesConst*
valueB :
џџџџџџџџџ*
dtype0*
_output_shapes
: 

metrics_1/mea/MeanMeanmetrics_1/mea/Abs$metrics_1/mea/Mean/reduction_indices*
T0*#
_output_shapes
:џџџџџџџџџ*
	keep_dims( *

Tidx0
]
metrics_1/mea/ConstConst*
valueB: *
dtype0*
_output_shapes
:

metrics_1/mea/SumSummetrics_1/mea/Meanmetrics_1/mea/Const*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
a
!metrics_1/mea/AssignAddVariableOpAssignAddVariableOptotal_2metrics_1/mea/Sum*
dtype0

metrics_1/mea/ReadVariableOpReadVariableOptotal_2"^metrics_1/mea/AssignAddVariableOp*
dtype0*
_output_shapes
: 
_
metrics_1/mea/SizeSizemetrics_1/mea/Mean*
T0*
out_type0*
_output_shapes
: 
n
metrics_1/mea/CastCastmetrics_1/mea/Size*

SrcT0*
Truncate( *
_output_shapes
: *

DstT0
d
#metrics_1/mea/AssignAddVariableOp_1AssignAddVariableOpcount_2metrics_1/mea/Cast*
dtype0

metrics_1/mea/ReadVariableOp_1ReadVariableOpcount_2$^metrics_1/mea/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ј
metrics_1/mea/ReadVariableOp_2ReadVariableOptotal_2"^metrics_1/mea/AssignAddVariableOp$^metrics_1/mea/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ў
$metrics_1/mea/truediv/ReadVariableOpReadVariableOpcount_2"^metrics_1/mea/AssignAddVariableOp$^metrics_1/mea/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

metrics_1/mea/truedivRealDivmetrics_1/mea/ReadVariableOp_2$metrics_1/mea/truediv/ReadVariableOp*
T0*
_output_shapes
: 
Z
metrics_1/mea/IdentityIdentitymetrics_1/mea/truediv*
_output_shapes
: *
T0
k
&metrics_1/msle/clip_by_value/Minimum/yConst*
valueB
 *  *
dtype0*
_output_shapes
: 
 
$metrics_1/msle/clip_by_value/MinimumMinimumactivation_8/Identity&metrics_1/msle/clip_by_value/Minimum/y*
T0*'
_output_shapes
:џџџџџџџџџ
c
metrics_1/msle/clip_by_value/yConst*
valueB
 *Пж3*
dtype0*
_output_shapes
: 

metrics_1/msle/clip_by_valueMaximum$metrics_1/msle/clip_by_value/Minimummetrics_1/msle/clip_by_value/y*
T0*'
_output_shapes
:џџџџџџџџџ
Y
metrics_1/msle/add/yConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 

metrics_1/msle/addAddmetrics_1/msle/clip_by_valuemetrics_1/msle/add/y*
T0*'
_output_shapes
:џџџџџџџџџ
_
metrics_1/msle/LogLogmetrics_1/msle/add*
T0*'
_output_shapes
:џџџџџџџџџ
m
(metrics_1/msle/clip_by_value_1/Minimum/yConst*
valueB
 *  *
dtype0*
_output_shapes
: 
Ћ
&metrics_1/msle/clip_by_value_1/MinimumMinimumactivation_8_target(metrics_1/msle/clip_by_value_1/Minimum/y*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
e
 metrics_1/msle/clip_by_value_1/yConst*
dtype0*
_output_shapes
: *
valueB
 *Пж3
Ў
metrics_1/msle/clip_by_value_1Maximum&metrics_1/msle/clip_by_value_1/Minimum metrics_1/msle/clip_by_value_1/y*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
[
metrics_1/msle/add_1/yConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 

metrics_1/msle/add_1Addmetrics_1/msle/clip_by_value_1metrics_1/msle/add_1/y*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
l
metrics_1/msle/Log_1Logmetrics_1/msle/add_1*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ*
T0
~
metrics_1/msle/subSubmetrics_1/msle/Logmetrics_1/msle/Log_1*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
n
metrics_1/msle/SquareSquaremetrics_1/msle/sub*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
p
%metrics_1/msle/Mean/reduction_indicesConst*
valueB :
џџџџџџџџџ*
dtype0*
_output_shapes
: 
Є
metrics_1/msle/MeanMeanmetrics_1/msle/Square%metrics_1/msle/Mean/reduction_indices*#
_output_shapes
:џџџџџџџџџ*
	keep_dims( *

Tidx0*
T0
^
metrics_1/msle/ConstConst*
valueB: *
dtype0*
_output_shapes
:

metrics_1/msle/SumSummetrics_1/msle/Meanmetrics_1/msle/Const*
	keep_dims( *

Tidx0*
T0*
_output_shapes
: 
c
"metrics_1/msle/AssignAddVariableOpAssignAddVariableOptotal_3metrics_1/msle/Sum*
dtype0

metrics_1/msle/ReadVariableOpReadVariableOptotal_3#^metrics_1/msle/AssignAddVariableOp*
dtype0*
_output_shapes
: 
a
metrics_1/msle/SizeSizemetrics_1/msle/Mean*
T0*
out_type0*
_output_shapes
: 
p
metrics_1/msle/CastCastmetrics_1/msle/Size*
Truncate( *
_output_shapes
: *

DstT0*

SrcT0
f
$metrics_1/msle/AssignAddVariableOp_1AssignAddVariableOpcount_3metrics_1/msle/Cast*
dtype0

metrics_1/msle/ReadVariableOp_1ReadVariableOpcount_3%^metrics_1/msle/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ћ
metrics_1/msle/ReadVariableOp_2ReadVariableOptotal_3#^metrics_1/msle/AssignAddVariableOp%^metrics_1/msle/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Б
%metrics_1/msle/truediv/ReadVariableOpReadVariableOpcount_3#^metrics_1/msle/AssignAddVariableOp%^metrics_1/msle/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

metrics_1/msle/truedivRealDivmetrics_1/msle/ReadVariableOp_2%metrics_1/msle/truediv/ReadVariableOp*
T0*
_output_shapes
: 
\
metrics_1/msle/IdentityIdentitymetrics_1/msle/truediv*
T0*
_output_shapes
: 

metrics_1/logcosh/subSubactivation_8/Identityactivation_8_target*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
\
metrics_1/logcosh/mul/xConst*
dtype0*
_output_shapes
: *
valueB
 *   Р

metrics_1/logcosh/mulMulmetrics_1/logcosh/mul/xmetrics_1/logcosh/sub*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
x
metrics_1/logcosh/SoftplusSoftplusmetrics_1/logcosh/mul*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ

metrics_1/logcosh/addAddmetrics_1/logcosh/submetrics_1/logcosh/Softplus*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
\
metrics_1/logcosh/Log/xConst*
valueB
 *   @*
dtype0*
_output_shapes
: 
V
metrics_1/logcosh/LogLogmetrics_1/logcosh/Log/x*
T0*
_output_shapes
: 

metrics_1/logcosh/sub_1Submetrics_1/logcosh/addmetrics_1/logcosh/Log*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
s
(metrics_1/logcosh/Mean/reduction_indicesConst*
dtype0*
_output_shapes
: *
valueB :
џџџџџџџџџ
Ќ
metrics_1/logcosh/MeanMeanmetrics_1/logcosh/sub_1(metrics_1/logcosh/Mean/reduction_indices*#
_output_shapes
:џџџџџџџџџ*
	keep_dims( *

Tidx0*
T0
a
metrics_1/logcosh/ConstConst*
valueB: *
dtype0*
_output_shapes
:

metrics_1/logcosh/SumSummetrics_1/logcosh/Meanmetrics_1/logcosh/Const*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
i
%metrics_1/logcosh/AssignAddVariableOpAssignAddVariableOptotal_4metrics_1/logcosh/Sum*
dtype0

 metrics_1/logcosh/ReadVariableOpReadVariableOptotal_4&^metrics_1/logcosh/AssignAddVariableOp*
dtype0*
_output_shapes
: 
g
metrics_1/logcosh/SizeSizemetrics_1/logcosh/Mean*
_output_shapes
: *
T0*
out_type0
v
metrics_1/logcosh/CastCastmetrics_1/logcosh/Size*

SrcT0*
Truncate( *
_output_shapes
: *

DstT0
l
'metrics_1/logcosh/AssignAddVariableOp_1AssignAddVariableOpcount_4metrics_1/logcosh/Cast*
dtype0

"metrics_1/logcosh/ReadVariableOp_1ReadVariableOpcount_4(^metrics_1/logcosh/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Д
"metrics_1/logcosh/ReadVariableOp_2ReadVariableOptotal_4&^metrics_1/logcosh/AssignAddVariableOp(^metrics_1/logcosh/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
К
(metrics_1/logcosh/truediv/ReadVariableOpReadVariableOpcount_4&^metrics_1/logcosh/AssignAddVariableOp(^metrics_1/logcosh/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

metrics_1/logcosh/truedivRealDiv"metrics_1/logcosh/ReadVariableOp_2(metrics_1/logcosh/truediv/ReadVariableOp*
T0*
_output_shapes
: 
b
metrics_1/logcosh/IdentityIdentitymetrics_1/logcosh/truediv*
T0*
_output_shapes
: 
e
metrics_1/corr/ConstConst*
dtype0*
_output_shapes
:*
valueB"       

metrics_1/corr/MeanMeanactivation_8_targetmetrics_1/corr/Const*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
g
metrics_1/corr/Const_1Const*
valueB"       *
dtype0*
_output_shapes
:

metrics_1/corr/Mean_1Meanactivation_8/Identitymetrics_1/corr/Const_1*
_output_shapes
: *
	keep_dims( *

Tidx0*
T0
~
metrics_1/corr/subSubactivation_8_targetmetrics_1/corr/Mean*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
{
metrics_1/corr/sub_1Subactivation_8/Identitymetrics_1/corr/Mean_1*
T0*'
_output_shapes
:џџџџџџџџџ
~
metrics_1/corr/mulMulmetrics_1/corr/submetrics_1/corr/sub_1*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
g
metrics_1/corr/Const_2Const*
valueB"       *
dtype0*
_output_shapes
:

metrics_1/corr/SumSummetrics_1/corr/mulmetrics_1/corr/Const_2*
_output_shapes
: *
	keep_dims( *

Tidx0*
T0
~
metrics_1/corr/mul_1Mulmetrics_1/corr/submetrics_1/corr/sub*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ
g
metrics_1/corr/Const_3Const*
dtype0*
_output_shapes
:*
valueB"       

metrics_1/corr/Sum_1Summetrics_1/corr/mul_1metrics_1/corr/Const_3*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
y
metrics_1/corr/mul_2Mulmetrics_1/corr/sub_1metrics_1/corr/sub_1*'
_output_shapes
:џџџџџџџџџ*
T0
g
metrics_1/corr/Const_4Const*
valueB"       *
dtype0*
_output_shapes
:

metrics_1/corr/Sum_2Summetrics_1/corr/mul_2metrics_1/corr/Const_4*
	keep_dims( *

Tidx0*
T0*
_output_shapes
: 
h
metrics_1/corr/mul_3Mulmetrics_1/corr/Sum_1metrics_1/corr/Sum_2*
T0*
_output_shapes
: 
[
metrics_1/corr/Const_5Const*
valueB
 *    *
dtype0*
_output_shapes
: 
[
metrics_1/corr/Const_6Const*
valueB
 *  *
dtype0*
_output_shapes
: 
~
$metrics_1/corr/clip_by_value/MinimumMinimummetrics_1/corr/mul_3metrics_1/corr/Const_6*
T0*
_output_shapes
: 

metrics_1/corr/clip_by_valueMaximum$metrics_1/corr/clip_by_value/Minimummetrics_1/corr/Const_5*
T0*
_output_shapes
: 
Z
metrics_1/corr/SqrtSqrtmetrics_1/corr/clip_by_value*
T0*
_output_shapes
: 
k
metrics_1/corr/truedivRealDivmetrics_1/corr/Summetrics_1/corr/Sqrt*
T0*
_output_shapes
: 
]
metrics_1/corr/Minimum/yConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
t
metrics_1/corr/MinimumMinimummetrics_1/corr/truedivmetrics_1/corr/Minimum/y*
_output_shapes
: *
T0
]
metrics_1/corr/Maximum/yConst*
dtype0*
_output_shapes
: *
valueB
 *  П
t
metrics_1/corr/MaximumMaximummetrics_1/corr/Minimummetrics_1/corr/Maximum/y*
_output_shapes
: *
T0
Y
metrics_1/corr/Const_7Const*
valueB *
dtype0*
_output_shapes
: 

metrics_1/corr/Sum_3Summetrics_1/corr/Maximummetrics_1/corr/Const_7*
	keep_dims( *

Tidx0*
T0*
_output_shapes
: 
e
"metrics_1/corr/AssignAddVariableOpAssignAddVariableOptotal_5metrics_1/corr/Sum_3*
dtype0

metrics_1/corr/ReadVariableOpReadVariableOptotal_5#^metrics_1/corr/AssignAddVariableOp*
dtype0*
_output_shapes
: 
U
metrics_1/corr/SizeConst*
value	B :*
dtype0*
_output_shapes
: 
p
metrics_1/corr/CastCastmetrics_1/corr/Size*

SrcT0*
Truncate( *
_output_shapes
: *

DstT0
f
$metrics_1/corr/AssignAddVariableOp_1AssignAddVariableOpcount_5metrics_1/corr/Cast*
dtype0

metrics_1/corr/ReadVariableOp_1ReadVariableOpcount_5%^metrics_1/corr/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Ћ
metrics_1/corr/ReadVariableOp_2ReadVariableOptotal_5#^metrics_1/corr/AssignAddVariableOp%^metrics_1/corr/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
Г
'metrics_1/corr/truediv_1/ReadVariableOpReadVariableOpcount_5#^metrics_1/corr/AssignAddVariableOp%^metrics_1/corr/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

metrics_1/corr/truediv_1RealDivmetrics_1/corr/ReadVariableOp_2'metrics_1/corr/truediv_1/ReadVariableOp*
T0*
_output_shapes
: 
^
metrics_1/corr/IdentityIdentitymetrics_1/corr/truediv_1*
_output_shapes
: *
T0

/loss_1/activation_8_loss/mean_squared_error/subSubactivation_8/Identityactivation_8_target*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ*
T0
Ј
2loss_1/activation_8_loss/mean_squared_error/SquareSquare/loss_1/activation_8_loss/mean_squared_error/sub*
T0*0
_output_shapes
:џџџџџџџџџџџџџџџџџџ

Bloss_1/activation_8_loss/mean_squared_error/Mean/reduction_indicesConst*
valueB :
џџџџџџџџџ*
dtype0*
_output_shapes
: 
ћ
0loss_1/activation_8_loss/mean_squared_error/MeanMean2loss_1/activation_8_loss/mean_squared_error/SquareBloss_1/activation_8_loss/mean_squared_error/Mean/reduction_indices*
	keep_dims( *

Tidx0*
T0*#
_output_shapes
:џџџџџџџџџ
С
=loss_1/activation_8_loss/mean_squared_error/weighted_loss/mulMulactivation_8_sample_weights0loss_1/activation_8_loss/mean_squared_error/Mean*#
_output_shapes
:џџџџџџџџџ*
T0

?loss_1/activation_8_loss/mean_squared_error/weighted_loss/ConstConst*
dtype0*
_output_shapes
:*
valueB: 

=loss_1/activation_8_loss/mean_squared_error/weighted_loss/SumSum=loss_1/activation_8_loss/mean_squared_error/weighted_loss/mul?loss_1/activation_8_loss/mean_squared_error/weighted_loss/Const*
T0*
_output_shapes
: *
	keep_dims( *

Tidx0
У
Kloss_1/activation_8_loss/mean_squared_error/weighted_loss/num_elements/SizeSize=loss_1/activation_8_loss/mean_squared_error/weighted_loss/mul*
_output_shapes
: *
T0*
out_type0
р
Kloss_1/activation_8_loss/mean_squared_error/weighted_loss/num_elements/CastCastKloss_1/activation_8_loss/mean_squared_error/weighted_loss/num_elements/Size*

SrcT0*
Truncate( *
_output_shapes
: *

DstT0
љ
Aloss_1/activation_8_loss/mean_squared_error/weighted_loss/truedivRealDiv=loss_1/activation_8_loss/mean_squared_error/weighted_loss/SumKloss_1/activation_8_loss/mean_squared_error/weighted_loss/num_elements/Cast*
T0*
_output_shapes
: 
Q
loss_1/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 


loss_1/mulMulloss_1/mul/xAloss_1/activation_8_loss/mean_squared_error/weighted_loss/truediv*
_output_shapes
: *
T0
K
Const_13Const*
valueB *
dtype0*
_output_shapes
: 
b
Mean_1Mean
loss_1/mulConst_13*
	keep_dims( *

Tidx0*
T0*
_output_shapes
: 
b
Placeholder_8Placeholder*
dtype0* 
_output_shapes
:
*
shape:

R
AssignVariableOp_8AssignVariableOpdense_5/kernelPlaceholder_8*
dtype0
v
ReadVariableOp_8ReadVariableOpdense_5/kernel^AssignVariableOp_8*
dtype0* 
_output_shapes
:

X
Placeholder_9Placeholder*
dtype0*
_output_shapes	
:*
shape:
P
AssignVariableOp_9AssignVariableOpdense_5/biasPlaceholder_9*
dtype0
o
ReadVariableOp_9ReadVariableOpdense_5/bias^AssignVariableOp_9*
dtype0*
_output_shapes	
:
c
Placeholder_10Placeholder*
dtype0* 
_output_shapes
:
*
shape:

T
AssignVariableOp_10AssignVariableOpdense_6/kernelPlaceholder_10*
dtype0
x
ReadVariableOp_10ReadVariableOpdense_6/kernel^AssignVariableOp_10*
dtype0* 
_output_shapes
:

Y
Placeholder_11Placeholder*
dtype0*
_output_shapes	
:*
shape:
R
AssignVariableOp_11AssignVariableOpdense_6/biasPlaceholder_11*
dtype0
q
ReadVariableOp_11ReadVariableOpdense_6/bias^AssignVariableOp_11*
dtype0*
_output_shapes	
:
c
Placeholder_12Placeholder*
dtype0* 
_output_shapes
:
*
shape:

T
AssignVariableOp_12AssignVariableOpdense_7/kernelPlaceholder_12*
dtype0
x
ReadVariableOp_12ReadVariableOpdense_7/kernel^AssignVariableOp_12*
dtype0* 
_output_shapes
:

Y
Placeholder_13Placeholder*
dtype0*
_output_shapes	
:*
shape:
R
AssignVariableOp_13AssignVariableOpdense_7/biasPlaceholder_13*
dtype0
q
ReadVariableOp_13ReadVariableOpdense_7/bias^AssignVariableOp_13*
dtype0*
_output_shapes	
:
a
Placeholder_14Placeholder*
dtype0*
_output_shapes
:	*
shape:	
T
AssignVariableOp_14AssignVariableOpdense_8/kernelPlaceholder_14*
dtype0
w
ReadVariableOp_14ReadVariableOpdense_8/kernel^AssignVariableOp_14*
dtype0*
_output_shapes
:	
W
Placeholder_15Placeholder*
shape:*
dtype0*
_output_shapes
:
R
AssignVariableOp_15AssignVariableOpdense_8/biasPlaceholder_15*
dtype0
p
ReadVariableOp_15ReadVariableOpdense_8/bias^AssignVariableOp_15*
dtype0*
_output_shapes
:
L
VarIsInitializedOp_15VarIsInitializedOpcount_3*
_output_shapes
: 
Q
VarIsInitializedOp_16VarIsInitializedOpdense_6/bias*
_output_shapes
: 
L
VarIsInitializedOp_17VarIsInitializedOpcount_4*
_output_shapes
: 
S
VarIsInitializedOp_18VarIsInitializedOpdense_7/kernel*
_output_shapes
: 
L
VarIsInitializedOp_19VarIsInitializedOptotal_5*
_output_shapes
: 
S
VarIsInitializedOp_20VarIsInitializedOpdense_6/kernel*
_output_shapes
: 
S
VarIsInitializedOp_21VarIsInitializedOpdense_5/kernel*
_output_shapes
: 
L
VarIsInitializedOp_22VarIsInitializedOptotal_4*
_output_shapes
: 
Q
VarIsInitializedOp_23VarIsInitializedOpdense_5/bias*
_output_shapes
: 
L
VarIsInitializedOp_24VarIsInitializedOpcount_5*
_output_shapes
: 
V
VarIsInitializedOp_25VarIsInitializedOpAdam_1/iterations*
_output_shapes
: 
L
VarIsInitializedOp_26VarIsInitializedOptotal_3*
_output_shapes
: 
S
VarIsInitializedOp_27VarIsInitializedOpdense_8/kernel*
_output_shapes
: 
Y
VarIsInitializedOp_28VarIsInitializedOpAdam_1/learning_rate*
_output_shapes
: 
R
VarIsInitializedOp_29VarIsInitializedOpAdam_1/beta_1*
_output_shapes
: 
Q
VarIsInitializedOp_30VarIsInitializedOpdense_7/bias*
_output_shapes
: 
R
VarIsInitializedOp_31VarIsInitializedOpAdam_1/beta_2*
_output_shapes
: 
Q
VarIsInitializedOp_32VarIsInitializedOpdense_8/bias*
_output_shapes
: 
Q
VarIsInitializedOp_33VarIsInitializedOpAdam_1/decay*
_output_shapes
: 
L
VarIsInitializedOp_34VarIsInitializedOptotal_1*
_output_shapes
: 
L
VarIsInitializedOp_35VarIsInitializedOptotal_2*
_output_shapes
: 
L
VarIsInitializedOp_36VarIsInitializedOpcount_1*
_output_shapes
: 
L
VarIsInitializedOp_37VarIsInitializedOpcount_2*
_output_shapes
: 
э
init_1NoOp^Adam_1/beta_1/Assign^Adam_1/beta_2/Assign^Adam_1/decay/Assign^Adam_1/iterations/Assign^Adam_1/learning_rate/Assign^count_1/Assign^count_2/Assign^count_3/Assign^count_4/Assign^count_5/Assign^dense_5/bias/Assign^dense_5/kernel/Assign^dense_6/bias/Assign^dense_6/kernel/Assign^dense_7/bias/Assign^dense_7/kernel/Assign^dense_8/bias/Assign^dense_8/kernel/Assign^total_1/Assign^total_2/Assign^total_3/Assign^total_4/Assign^total_5/Assign
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
dtype0*
_output_shapes
: *
shape: 
e

save/ConstPlaceholderWithDefaultsave/filename*
dtype0*
_output_shapes
: *
shape: 

save/StringJoin/inputs_1Const*<
value3B1 B+_temp_05ae730245a44cc6bd9f07dc9a8f75d5/part*
dtype0*
_output_shapes
: 
u
save/StringJoin
StringJoin
save/Constsave/StringJoin/inputs_1*
N*
_output_shapes
: *
	separator 
Q
save/num_shardsConst*
value	B :*
dtype0*
_output_shapes
: 
k
save/ShardedFilename/shardConst"/device:CPU:0*
value	B : *
dtype0*
_output_shapes
: 

save/ShardedFilenameShardedFilenamesave/StringJoinsave/ShardedFilename/shardsave/num_shards"/device:CPU:0*
_output_shapes
: 
щ
save/SaveV2/tensor_namesConst"/device:CPU:0*
valueB&BAdam/beta_1BAdam/beta_2B
Adam/decayBAdam/iterationsBAdam/learning_rateBAdam_1/beta_1BAdam_1/beta_2BAdam_1/decayBAdam_1/iterationsBAdam_1/learning_rateBcountBcount_1Bcount_2Bcount_3Bcount_4Bcount_5Bdense_1/biasBdense_1/kernelBdense_2/biasBdense_2/kernelBdense_3/biasBdense_3/kernelBdense_4/biasBdense_4/kernelBdense_5/biasBdense_5/kernelBdense_6/biasBdense_6/kernelBdense_7/biasBdense_7/kernelBdense_8/biasBdense_8/kernelBtotalBtotal_1Btotal_2Btotal_3Btotal_4Btotal_5*
dtype0*
_output_shapes
:&
О
save/SaveV2/shape_and_slicesConst"/device:CPU:0*_
valueVBT&B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:&

save/SaveV2SaveV2save/ShardedFilenamesave/SaveV2/tensor_namessave/SaveV2/shape_and_slicesAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp#Adam/iterations/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOp!Adam_1/beta_1/Read/ReadVariableOp!Adam_1/beta_2/Read/ReadVariableOp Adam_1/decay/Read/ReadVariableOp%Adam_1/iterations/Read/ReadVariableOp(Adam_1/learning_rate/Read/ReadVariableOpcount/Read/ReadVariableOpcount_1/Read/ReadVariableOpcount_2/Read/ReadVariableOpcount_3/Read/ReadVariableOpcount_4/Read/ReadVariableOpcount_5/Read/ReadVariableOp dense_1/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOp"dense_2/kernel/Read/ReadVariableOp dense_3/bias/Read/ReadVariableOp"dense_3/kernel/Read/ReadVariableOp dense_4/bias/Read/ReadVariableOp"dense_4/kernel/Read/ReadVariableOp dense_5/bias/Read/ReadVariableOp"dense_5/kernel/Read/ReadVariableOp dense_6/bias/Read/ReadVariableOp"dense_6/kernel/Read/ReadVariableOp dense_7/bias/Read/ReadVariableOp"dense_7/kernel/Read/ReadVariableOp dense_8/bias/Read/ReadVariableOp"dense_8/kernel/Read/ReadVariableOptotal/Read/ReadVariableOptotal_1/Read/ReadVariableOptotal_2/Read/ReadVariableOptotal_3/Read/ReadVariableOptotal_4/Read/ReadVariableOptotal_5/Read/ReadVariableOp"/device:CPU:0*4
dtypes*
(2&		
 
save/control_dependencyIdentitysave/ShardedFilename^save/SaveV2"/device:CPU:0*
T0*'
_class
loc:@save/ShardedFilename*
_output_shapes
: 
Ќ
+save/MergeV2Checkpoints/checkpoint_prefixesPacksave/ShardedFilename^save/control_dependency"/device:CPU:0*
T0*

axis *
N*
_output_shapes
:

save/MergeV2CheckpointsMergeV2Checkpoints+save/MergeV2Checkpoints/checkpoint_prefixes
save/Const"/device:CPU:0*
delete_old_dirs(

save/IdentityIdentity
save/Const^save/MergeV2Checkpoints^save/control_dependency"/device:CPU:0*
T0*
_output_shapes
: 
ь
save/RestoreV2/tensor_namesConst"/device:CPU:0*
valueB&BAdam/beta_1BAdam/beta_2B
Adam/decayBAdam/iterationsBAdam/learning_rateBAdam_1/beta_1BAdam_1/beta_2BAdam_1/decayBAdam_1/iterationsBAdam_1/learning_rateBcountBcount_1Bcount_2Bcount_3Bcount_4Bcount_5Bdense_1/biasBdense_1/kernelBdense_2/biasBdense_2/kernelBdense_3/biasBdense_3/kernelBdense_4/biasBdense_4/kernelBdense_5/biasBdense_5/kernelBdense_6/biasBdense_6/kernelBdense_7/biasBdense_7/kernelBdense_8/biasBdense_8/kernelBtotalBtotal_1Btotal_2Btotal_3Btotal_4Btotal_5*
dtype0*
_output_shapes
:&
С
save/RestoreV2/shape_and_slicesConst"/device:CPU:0*_
valueVBT&B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:&
л
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*4
dtypes*
(2&		*Ў
_output_shapes
::::::::::::::::::::::::::::::::::::::
N
save/Identity_1Identitysave/RestoreV2*
T0*
_output_shapes
:
T
save/AssignVariableOpAssignVariableOpAdam/beta_1save/Identity_1*
dtype0
P
save/Identity_2Identitysave/RestoreV2:1*
T0*
_output_shapes
:
V
save/AssignVariableOp_1AssignVariableOpAdam/beta_2save/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:2*
T0*
_output_shapes
:
U
save/AssignVariableOp_2AssignVariableOp
Adam/decaysave/Identity_3*
dtype0
P
save/Identity_4Identitysave/RestoreV2:3*
_output_shapes
:*
T0	
Z
save/AssignVariableOp_3AssignVariableOpAdam/iterationssave/Identity_4*
dtype0	
P
save/Identity_5Identitysave/RestoreV2:4*
T0*
_output_shapes
:
]
save/AssignVariableOp_4AssignVariableOpAdam/learning_ratesave/Identity_5*
dtype0
P
save/Identity_6Identitysave/RestoreV2:5*
T0*
_output_shapes
:
X
save/AssignVariableOp_5AssignVariableOpAdam_1/beta_1save/Identity_6*
dtype0
P
save/Identity_7Identitysave/RestoreV2:6*
T0*
_output_shapes
:
X
save/AssignVariableOp_6AssignVariableOpAdam_1/beta_2save/Identity_7*
dtype0
P
save/Identity_8Identitysave/RestoreV2:7*
T0*
_output_shapes
:
W
save/AssignVariableOp_7AssignVariableOpAdam_1/decaysave/Identity_8*
dtype0
P
save/Identity_9Identitysave/RestoreV2:8*
T0	*
_output_shapes
:
\
save/AssignVariableOp_8AssignVariableOpAdam_1/iterationssave/Identity_9*
dtype0	
Q
save/Identity_10Identitysave/RestoreV2:9*
_output_shapes
:*
T0
`
save/AssignVariableOp_9AssignVariableOpAdam_1/learning_ratesave/Identity_10*
dtype0
R
save/Identity_11Identitysave/RestoreV2:10*
_output_shapes
:*
T0
R
save/AssignVariableOp_10AssignVariableOpcountsave/Identity_11*
dtype0
R
save/Identity_12Identitysave/RestoreV2:11*
T0*
_output_shapes
:
T
save/AssignVariableOp_11AssignVariableOpcount_1save/Identity_12*
dtype0
R
save/Identity_13Identitysave/RestoreV2:12*
T0*
_output_shapes
:
T
save/AssignVariableOp_12AssignVariableOpcount_2save/Identity_13*
dtype0
R
save/Identity_14Identitysave/RestoreV2:13*
T0*
_output_shapes
:
T
save/AssignVariableOp_13AssignVariableOpcount_3save/Identity_14*
dtype0
R
save/Identity_15Identitysave/RestoreV2:14*
_output_shapes
:*
T0
T
save/AssignVariableOp_14AssignVariableOpcount_4save/Identity_15*
dtype0
R
save/Identity_16Identitysave/RestoreV2:15*
T0*
_output_shapes
:
T
save/AssignVariableOp_15AssignVariableOpcount_5save/Identity_16*
dtype0
R
save/Identity_17Identitysave/RestoreV2:16*
T0*
_output_shapes
:
Y
save/AssignVariableOp_16AssignVariableOpdense_1/biassave/Identity_17*
dtype0
R
save/Identity_18Identitysave/RestoreV2:17*
T0*
_output_shapes
:
[
save/AssignVariableOp_17AssignVariableOpdense_1/kernelsave/Identity_18*
dtype0
R
save/Identity_19Identitysave/RestoreV2:18*
T0*
_output_shapes
:
Y
save/AssignVariableOp_18AssignVariableOpdense_2/biassave/Identity_19*
dtype0
R
save/Identity_20Identitysave/RestoreV2:19*
_output_shapes
:*
T0
[
save/AssignVariableOp_19AssignVariableOpdense_2/kernelsave/Identity_20*
dtype0
R
save/Identity_21Identitysave/RestoreV2:20*
T0*
_output_shapes
:
Y
save/AssignVariableOp_20AssignVariableOpdense_3/biassave/Identity_21*
dtype0
R
save/Identity_22Identitysave/RestoreV2:21*
T0*
_output_shapes
:
[
save/AssignVariableOp_21AssignVariableOpdense_3/kernelsave/Identity_22*
dtype0
R
save/Identity_23Identitysave/RestoreV2:22*
T0*
_output_shapes
:
Y
save/AssignVariableOp_22AssignVariableOpdense_4/biassave/Identity_23*
dtype0
R
save/Identity_24Identitysave/RestoreV2:23*
T0*
_output_shapes
:
[
save/AssignVariableOp_23AssignVariableOpdense_4/kernelsave/Identity_24*
dtype0
R
save/Identity_25Identitysave/RestoreV2:24*
T0*
_output_shapes
:
Y
save/AssignVariableOp_24AssignVariableOpdense_5/biassave/Identity_25*
dtype0
R
save/Identity_26Identitysave/RestoreV2:25*
T0*
_output_shapes
:
[
save/AssignVariableOp_25AssignVariableOpdense_5/kernelsave/Identity_26*
dtype0
R
save/Identity_27Identitysave/RestoreV2:26*
_output_shapes
:*
T0
Y
save/AssignVariableOp_26AssignVariableOpdense_6/biassave/Identity_27*
dtype0
R
save/Identity_28Identitysave/RestoreV2:27*
T0*
_output_shapes
:
[
save/AssignVariableOp_27AssignVariableOpdense_6/kernelsave/Identity_28*
dtype0
R
save/Identity_29Identitysave/RestoreV2:28*
T0*
_output_shapes
:
Y
save/AssignVariableOp_28AssignVariableOpdense_7/biassave/Identity_29*
dtype0
R
save/Identity_30Identitysave/RestoreV2:29*
_output_shapes
:*
T0
[
save/AssignVariableOp_29AssignVariableOpdense_7/kernelsave/Identity_30*
dtype0
R
save/Identity_31Identitysave/RestoreV2:30*
T0*
_output_shapes
:
Y
save/AssignVariableOp_30AssignVariableOpdense_8/biassave/Identity_31*
dtype0
R
save/Identity_32Identitysave/RestoreV2:31*
_output_shapes
:*
T0
[
save/AssignVariableOp_31AssignVariableOpdense_8/kernelsave/Identity_32*
dtype0
R
save/Identity_33Identitysave/RestoreV2:32*
T0*
_output_shapes
:
R
save/AssignVariableOp_32AssignVariableOptotalsave/Identity_33*
dtype0
R
save/Identity_34Identitysave/RestoreV2:33*
T0*
_output_shapes
:
T
save/AssignVariableOp_33AssignVariableOptotal_1save/Identity_34*
dtype0
R
save/Identity_35Identitysave/RestoreV2:34*
T0*
_output_shapes
:
T
save/AssignVariableOp_34AssignVariableOptotal_2save/Identity_35*
dtype0
R
save/Identity_36Identitysave/RestoreV2:35*
T0*
_output_shapes
:
T
save/AssignVariableOp_35AssignVariableOptotal_3save/Identity_36*
dtype0
R
save/Identity_37Identitysave/RestoreV2:36*
T0*
_output_shapes
:
T
save/AssignVariableOp_36AssignVariableOptotal_4save/Identity_37*
dtype0
R
save/Identity_38Identitysave/RestoreV2:37*
T0*
_output_shapes
:
T
save/AssignVariableOp_37AssignVariableOptotal_5save/Identity_38*
dtype0

save/restore_shardNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_10^save/AssignVariableOp_11^save/AssignVariableOp_12^save/AssignVariableOp_13^save/AssignVariableOp_14^save/AssignVariableOp_15^save/AssignVariableOp_16^save/AssignVariableOp_17^save/AssignVariableOp_18^save/AssignVariableOp_19^save/AssignVariableOp_2^save/AssignVariableOp_20^save/AssignVariableOp_21^save/AssignVariableOp_22^save/AssignVariableOp_23^save/AssignVariableOp_24^save/AssignVariableOp_25^save/AssignVariableOp_26^save/AssignVariableOp_27^save/AssignVariableOp_28^save/AssignVariableOp_29^save/AssignVariableOp_3^save/AssignVariableOp_30^save/AssignVariableOp_31^save/AssignVariableOp_32^save/AssignVariableOp_33^save/AssignVariableOp_34^save/AssignVariableOp_35^save/AssignVariableOp_36^save/AssignVariableOp_37^save/AssignVariableOp_4^save/AssignVariableOp_5^save/AssignVariableOp_6^save/AssignVariableOp_7^save/AssignVariableOp_8^save/AssignVariableOp_9
-
save/restore_allNoOp^save/restore_shard"&<
save/Const:0save/Identity:0save/restore_all (5 @F8"ЖC
cond_contextЅCЂC

dropout_1/cond/cond_textdropout_1/cond/pred_id:0dropout_1/cond/switch_t:0 *М
dense_1/BiasAdd:0
dropout_1/cond/dropout/Cast:0
%dropout_1/cond/dropout/GreaterEqual:0
%dropout_1/cond/dropout/Shape/Switch:1
dropout_1/cond/dropout/Shape:0
dropout_1/cond/dropout/mul:0
dropout_1/cond/dropout/mul_1:0
5dropout_1/cond/dropout/random_uniform/RandomUniform:0
+dropout_1/cond/dropout/random_uniform/max:0
+dropout_1/cond/dropout/random_uniform/min:0
+dropout_1/cond/dropout/random_uniform/mul:0
+dropout_1/cond/dropout/random_uniform/sub:0
'dropout_1/cond/dropout/random_uniform:0
dropout_1/cond/dropout/rate:0
dropout_1/cond/dropout/sub/x:0
dropout_1/cond/dropout/sub:0
"dropout_1/cond/dropout/truediv/x:0
 dropout_1/cond/dropout/truediv:0
dropout_1/cond/pred_id:0
dropout_1/cond/switch_t:04
dropout_1/cond/pred_id:0dropout_1/cond/pred_id:0:
dense_1/BiasAdd:0%dropout_1/cond/dropout/Shape/Switch:1
И
dropout_1/cond/cond_text_1dropout_1/cond/pred_id:0dropout_1/cond/switch_f:0*ф
dense_1/BiasAdd:0
dropout_1/cond/Switch_1:0
dropout_1/cond/Switch_1:1
dropout_1/cond/pred_id:0
dropout_1/cond/switch_f:04
dropout_1/cond/pred_id:0dropout_1/cond/pred_id:0.
dense_1/BiasAdd:0dropout_1/cond/Switch_1:0

dropout_2/cond/cond_textdropout_2/cond/pred_id:0dropout_2/cond/switch_t:0 *М
dense_2/BiasAdd:0
dropout_2/cond/dropout/Cast:0
%dropout_2/cond/dropout/GreaterEqual:0
%dropout_2/cond/dropout/Shape/Switch:1
dropout_2/cond/dropout/Shape:0
dropout_2/cond/dropout/mul:0
dropout_2/cond/dropout/mul_1:0
5dropout_2/cond/dropout/random_uniform/RandomUniform:0
+dropout_2/cond/dropout/random_uniform/max:0
+dropout_2/cond/dropout/random_uniform/min:0
+dropout_2/cond/dropout/random_uniform/mul:0
+dropout_2/cond/dropout/random_uniform/sub:0
'dropout_2/cond/dropout/random_uniform:0
dropout_2/cond/dropout/rate:0
dropout_2/cond/dropout/sub/x:0
dropout_2/cond/dropout/sub:0
"dropout_2/cond/dropout/truediv/x:0
 dropout_2/cond/dropout/truediv:0
dropout_2/cond/pred_id:0
dropout_2/cond/switch_t:0:
dense_2/BiasAdd:0%dropout_2/cond/dropout/Shape/Switch:14
dropout_2/cond/pred_id:0dropout_2/cond/pred_id:0
И
dropout_2/cond/cond_text_1dropout_2/cond/pred_id:0dropout_2/cond/switch_f:0*ф
dense_2/BiasAdd:0
dropout_2/cond/Switch_1:0
dropout_2/cond/Switch_1:1
dropout_2/cond/pred_id:0
dropout_2/cond/switch_f:0.
dense_2/BiasAdd:0dropout_2/cond/Switch_1:04
dropout_2/cond/pred_id:0dropout_2/cond/pred_id:0

dropout_3/cond/cond_textdropout_3/cond/pred_id:0dropout_3/cond/switch_t:0 *М
dense_3/BiasAdd:0
dropout_3/cond/dropout/Cast:0
%dropout_3/cond/dropout/GreaterEqual:0
%dropout_3/cond/dropout/Shape/Switch:1
dropout_3/cond/dropout/Shape:0
dropout_3/cond/dropout/mul:0
dropout_3/cond/dropout/mul_1:0
5dropout_3/cond/dropout/random_uniform/RandomUniform:0
+dropout_3/cond/dropout/random_uniform/max:0
+dropout_3/cond/dropout/random_uniform/min:0
+dropout_3/cond/dropout/random_uniform/mul:0
+dropout_3/cond/dropout/random_uniform/sub:0
'dropout_3/cond/dropout/random_uniform:0
dropout_3/cond/dropout/rate:0
dropout_3/cond/dropout/sub/x:0
dropout_3/cond/dropout/sub:0
"dropout_3/cond/dropout/truediv/x:0
 dropout_3/cond/dropout/truediv:0
dropout_3/cond/pred_id:0
dropout_3/cond/switch_t:04
dropout_3/cond/pred_id:0dropout_3/cond/pred_id:0:
dense_3/BiasAdd:0%dropout_3/cond/dropout/Shape/Switch:1
И
dropout_3/cond/cond_text_1dropout_3/cond/pred_id:0dropout_3/cond/switch_f:0*ф
dense_3/BiasAdd:0
dropout_3/cond/Switch_1:0
dropout_3/cond/Switch_1:1
dropout_3/cond/pred_id:0
dropout_3/cond/switch_f:04
dropout_3/cond/pred_id:0dropout_3/cond/pred_id:0.
dense_3/BiasAdd:0dropout_3/cond/Switch_1:0

dropout_4/cond/cond_textdropout_4/cond/pred_id:0dropout_4/cond/switch_t:0 *М
dropout_4/cond/dropout/Cast:0
%dropout_4/cond/dropout/GreaterEqual:0
%dropout_4/cond/dropout/Shape/Switch:1
dropout_4/cond/dropout/Shape:0
dropout_4/cond/dropout/mul:0
dropout_4/cond/dropout/mul_1:0
5dropout_4/cond/dropout/random_uniform/RandomUniform:0
+dropout_4/cond/dropout/random_uniform/max:0
+dropout_4/cond/dropout/random_uniform/min:0
+dropout_4/cond/dropout/random_uniform/mul:0
+dropout_4/cond/dropout/random_uniform/sub:0
'dropout_4/cond/dropout/random_uniform:0
dropout_4/cond/dropout/rate:0
dropout_4/cond/dropout/sub/x:0
dropout_4/cond/dropout/sub:0
"dropout_4/cond/dropout/truediv/x:0
 dropout_4/cond/dropout/truediv:0
dropout_4/cond/pred_id:0
dropout_4/cond/switch_t:0
dropout_4_input:0:
dropout_4_input:0%dropout_4/cond/dropout/Shape/Switch:14
dropout_4/cond/pred_id:0dropout_4/cond/pred_id:0
И
dropout_4/cond/cond_text_1dropout_4/cond/pred_id:0dropout_4/cond/switch_f:0*ф
dropout_4/cond/Switch_1:0
dropout_4/cond/Switch_1:1
dropout_4/cond/pred_id:0
dropout_4/cond/switch_f:0
dropout_4_input:04
dropout_4/cond/pred_id:0dropout_4/cond/pred_id:0.
dropout_4_input:0dropout_4/cond/Switch_1:0

dropout_5/cond/cond_textdropout_5/cond/pred_id:0dropout_5/cond/switch_t:0 *М
dense_5/BiasAdd:0
dropout_5/cond/dropout/Cast:0
%dropout_5/cond/dropout/GreaterEqual:0
%dropout_5/cond/dropout/Shape/Switch:1
dropout_5/cond/dropout/Shape:0
dropout_5/cond/dropout/mul:0
dropout_5/cond/dropout/mul_1:0
5dropout_5/cond/dropout/random_uniform/RandomUniform:0
+dropout_5/cond/dropout/random_uniform/max:0
+dropout_5/cond/dropout/random_uniform/min:0
+dropout_5/cond/dropout/random_uniform/mul:0
+dropout_5/cond/dropout/random_uniform/sub:0
'dropout_5/cond/dropout/random_uniform:0
dropout_5/cond/dropout/rate:0
dropout_5/cond/dropout/sub/x:0
dropout_5/cond/dropout/sub:0
"dropout_5/cond/dropout/truediv/x:0
 dropout_5/cond/dropout/truediv:0
dropout_5/cond/pred_id:0
dropout_5/cond/switch_t:04
dropout_5/cond/pred_id:0dropout_5/cond/pred_id:0:
dense_5/BiasAdd:0%dropout_5/cond/dropout/Shape/Switch:1
И
dropout_5/cond/cond_text_1dropout_5/cond/pred_id:0dropout_5/cond/switch_f:0*ф
dense_5/BiasAdd:0
dropout_5/cond/Switch_1:0
dropout_5/cond/Switch_1:1
dropout_5/cond/pred_id:0
dropout_5/cond/switch_f:04
dropout_5/cond/pred_id:0dropout_5/cond/pred_id:0.
dense_5/BiasAdd:0dropout_5/cond/Switch_1:0

dropout_6/cond/cond_textdropout_6/cond/pred_id:0dropout_6/cond/switch_t:0 *М
dense_6/BiasAdd:0
dropout_6/cond/dropout/Cast:0
%dropout_6/cond/dropout/GreaterEqual:0
%dropout_6/cond/dropout/Shape/Switch:1
dropout_6/cond/dropout/Shape:0
dropout_6/cond/dropout/mul:0
dropout_6/cond/dropout/mul_1:0
5dropout_6/cond/dropout/random_uniform/RandomUniform:0
+dropout_6/cond/dropout/random_uniform/max:0
+dropout_6/cond/dropout/random_uniform/min:0
+dropout_6/cond/dropout/random_uniform/mul:0
+dropout_6/cond/dropout/random_uniform/sub:0
'dropout_6/cond/dropout/random_uniform:0
dropout_6/cond/dropout/rate:0
dropout_6/cond/dropout/sub/x:0
dropout_6/cond/dropout/sub:0
"dropout_6/cond/dropout/truediv/x:0
 dropout_6/cond/dropout/truediv:0
dropout_6/cond/pred_id:0
dropout_6/cond/switch_t:0:
dense_6/BiasAdd:0%dropout_6/cond/dropout/Shape/Switch:14
dropout_6/cond/pred_id:0dropout_6/cond/pred_id:0
И
dropout_6/cond/cond_text_1dropout_6/cond/pred_id:0dropout_6/cond/switch_f:0*ф
dense_6/BiasAdd:0
dropout_6/cond/Switch_1:0
dropout_6/cond/Switch_1:1
dropout_6/cond/pred_id:0
dropout_6/cond/switch_f:0.
dense_6/BiasAdd:0dropout_6/cond/Switch_1:04
dropout_6/cond/pred_id:0dropout_6/cond/pred_id:0

dropout_7/cond/cond_textdropout_7/cond/pred_id:0dropout_7/cond/switch_t:0 *М
dense_7/BiasAdd:0
dropout_7/cond/dropout/Cast:0
%dropout_7/cond/dropout/GreaterEqual:0
%dropout_7/cond/dropout/Shape/Switch:1
dropout_7/cond/dropout/Shape:0
dropout_7/cond/dropout/mul:0
dropout_7/cond/dropout/mul_1:0
5dropout_7/cond/dropout/random_uniform/RandomUniform:0
+dropout_7/cond/dropout/random_uniform/max:0
+dropout_7/cond/dropout/random_uniform/min:0
+dropout_7/cond/dropout/random_uniform/mul:0
+dropout_7/cond/dropout/random_uniform/sub:0
'dropout_7/cond/dropout/random_uniform:0
dropout_7/cond/dropout/rate:0
dropout_7/cond/dropout/sub/x:0
dropout_7/cond/dropout/sub:0
"dropout_7/cond/dropout/truediv/x:0
 dropout_7/cond/dropout/truediv:0
dropout_7/cond/pred_id:0
dropout_7/cond/switch_t:04
dropout_7/cond/pred_id:0dropout_7/cond/pred_id:0:
dense_7/BiasAdd:0%dropout_7/cond/dropout/Shape/Switch:1
И
dropout_7/cond/cond_text_1dropout_7/cond/pred_id:0dropout_7/cond/switch_f:0*ф
dense_7/BiasAdd:0
dropout_7/cond/Switch_1:0
dropout_7/cond/Switch_1:1
dropout_7/cond/pred_id:0
dropout_7/cond/switch_f:04
dropout_7/cond/pred_id:0dropout_7/cond/pred_id:0.
dense_7/BiasAdd:0dropout_7/cond/Switch_1:0"
	variablesџ
m
dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2dense_1/random_uniform:08
^
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2dense_1/Const:08
m
dense_2/kernel:0dense_2/kernel/Assign$dense_2/kernel/Read/ReadVariableOp:0(2dense_2/random_uniform:08
^
dense_2/bias:0dense_2/bias/Assign"dense_2/bias/Read/ReadVariableOp:0(2dense_2/Const:08
m
dense_3/kernel:0dense_3/kernel/Assign$dense_3/kernel/Read/ReadVariableOp:0(2dense_3/random_uniform:08
^
dense_3/bias:0dense_3/bias/Assign"dense_3/bias/Read/ReadVariableOp:0(2dense_3/Const:08
m
dense_4/kernel:0dense_4/kernel/Assign$dense_4/kernel/Read/ReadVariableOp:0(2dense_4/random_uniform:08
^
dense_4/bias:0dense_4/bias/Assign"dense_4/bias/Read/ReadVariableOp:0(2dense_4/Const:08

Adam/iterations:0Adam/iterations/Assign%Adam/iterations/Read/ReadVariableOp:0(2+Adam/iterations/Initializer/initial_value:08

Adam/learning_rate:0Adam/learning_rate/Assign(Adam/learning_rate/Read/ReadVariableOp:0(2.Adam/learning_rate/Initializer/initial_value:08
s
Adam/beta_1:0Adam/beta_1/Assign!Adam/beta_1/Read/ReadVariableOp:0(2'Adam/beta_1/Initializer/initial_value:08
s
Adam/beta_2:0Adam/beta_2/Assign!Adam/beta_2/Read/ReadVariableOp:0(2'Adam/beta_2/Initializer/initial_value:08
o
Adam/decay:0Adam/decay/Assign Adam/decay/Read/ReadVariableOp:0(2&Adam/decay/Initializer/initial_value:08
A
total:0total/Assigntotal/Read/ReadVariableOp:0(2Const:08
C
count:0count/Assigncount/Read/ReadVariableOp:0(2	Const_1:08
m
dense_5/kernel:0dense_5/kernel/Assign$dense_5/kernel/Read/ReadVariableOp:0(2dense_5/random_uniform:08
^
dense_5/bias:0dense_5/bias/Assign"dense_5/bias/Read/ReadVariableOp:0(2dense_5/Const:08
m
dense_6/kernel:0dense_6/kernel/Assign$dense_6/kernel/Read/ReadVariableOp:0(2dense_6/random_uniform:08
^
dense_6/bias:0dense_6/bias/Assign"dense_6/bias/Read/ReadVariableOp:0(2dense_6/Const:08
m
dense_7/kernel:0dense_7/kernel/Assign$dense_7/kernel/Read/ReadVariableOp:0(2dense_7/random_uniform:08
^
dense_7/bias:0dense_7/bias/Assign"dense_7/bias/Read/ReadVariableOp:0(2dense_7/Const:08
m
dense_8/kernel:0dense_8/kernel/Assign$dense_8/kernel/Read/ReadVariableOp:0(2dense_8/random_uniform:08
^
dense_8/bias:0dense_8/bias/Assign"dense_8/bias/Read/ReadVariableOp:0(2dense_8/Const:08

Adam_1/iterations:0Adam_1/iterations/Assign'Adam_1/iterations/Read/ReadVariableOp:0(2-Adam_1/iterations/Initializer/initial_value:08

Adam_1/learning_rate:0Adam_1/learning_rate/Assign*Adam_1/learning_rate/Read/ReadVariableOp:0(20Adam_1/learning_rate/Initializer/initial_value:08
{
Adam_1/beta_1:0Adam_1/beta_1/Assign#Adam_1/beta_1/Read/ReadVariableOp:0(2)Adam_1/beta_1/Initializer/initial_value:08
{
Adam_1/beta_2:0Adam_1/beta_2/Assign#Adam_1/beta_2/Read/ReadVariableOp:0(2)Adam_1/beta_2/Initializer/initial_value:08
w
Adam_1/decay:0Adam_1/decay/Assign"Adam_1/decay/Read/ReadVariableOp:0(2(Adam_1/decay/Initializer/initial_value:08
I
	total_1:0total_1/Assigntotal_1/Read/ReadVariableOp:0(2	Const_3:08
I
	count_1:0count_1/Assigncount_1/Read/ReadVariableOp:0(2	Const_4:08
I
	total_2:0total_2/Assigntotal_2/Read/ReadVariableOp:0(2	Const_5:08
I
	count_2:0count_2/Assigncount_2/Read/ReadVariableOp:0(2	Const_6:08
I
	total_3:0total_3/Assigntotal_3/Read/ReadVariableOp:0(2	Const_7:08
I
	count_3:0count_3/Assigncount_3/Read/ReadVariableOp:0(2	Const_8:08
I
	total_4:0total_4/Assigntotal_4/Read/ReadVariableOp:0(2	Const_9:08
J
	count_4:0count_4/Assigncount_4/Read/ReadVariableOp:0(2
Const_10:08
J
	total_5:0total_5/Assigntotal_5/Read/ReadVariableOp:0(2
Const_11:08
J
	count_5:0count_5/Assigncount_5/Read/ReadVariableOp:0(2
Const_12:08"
trainable_variablesџ
m
dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2dense_1/random_uniform:08
^
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2dense_1/Const:08
m
dense_2/kernel:0dense_2/kernel/Assign$dense_2/kernel/Read/ReadVariableOp:0(2dense_2/random_uniform:08
^
dense_2/bias:0dense_2/bias/Assign"dense_2/bias/Read/ReadVariableOp:0(2dense_2/Const:08
m
dense_3/kernel:0dense_3/kernel/Assign$dense_3/kernel/Read/ReadVariableOp:0(2dense_3/random_uniform:08
^
dense_3/bias:0dense_3/bias/Assign"dense_3/bias/Read/ReadVariableOp:0(2dense_3/Const:08
m
dense_4/kernel:0dense_4/kernel/Assign$dense_4/kernel/Read/ReadVariableOp:0(2dense_4/random_uniform:08
^
dense_4/bias:0dense_4/bias/Assign"dense_4/bias/Read/ReadVariableOp:0(2dense_4/Const:08

Adam/iterations:0Adam/iterations/Assign%Adam/iterations/Read/ReadVariableOp:0(2+Adam/iterations/Initializer/initial_value:08

Adam/learning_rate:0Adam/learning_rate/Assign(Adam/learning_rate/Read/ReadVariableOp:0(2.Adam/learning_rate/Initializer/initial_value:08
s
Adam/beta_1:0Adam/beta_1/Assign!Adam/beta_1/Read/ReadVariableOp:0(2'Adam/beta_1/Initializer/initial_value:08
s
Adam/beta_2:0Adam/beta_2/Assign!Adam/beta_2/Read/ReadVariableOp:0(2'Adam/beta_2/Initializer/initial_value:08
o
Adam/decay:0Adam/decay/Assign Adam/decay/Read/ReadVariableOp:0(2&Adam/decay/Initializer/initial_value:08
A
total:0total/Assigntotal/Read/ReadVariableOp:0(2Const:08
C
count:0count/Assigncount/Read/ReadVariableOp:0(2	Const_1:08
m
dense_5/kernel:0dense_5/kernel/Assign$dense_5/kernel/Read/ReadVariableOp:0(2dense_5/random_uniform:08
^
dense_5/bias:0dense_5/bias/Assign"dense_5/bias/Read/ReadVariableOp:0(2dense_5/Const:08
m
dense_6/kernel:0dense_6/kernel/Assign$dense_6/kernel/Read/ReadVariableOp:0(2dense_6/random_uniform:08
^
dense_6/bias:0dense_6/bias/Assign"dense_6/bias/Read/ReadVariableOp:0(2dense_6/Const:08
m
dense_7/kernel:0dense_7/kernel/Assign$dense_7/kernel/Read/ReadVariableOp:0(2dense_7/random_uniform:08
^
dense_7/bias:0dense_7/bias/Assign"dense_7/bias/Read/ReadVariableOp:0(2dense_7/Const:08
m
dense_8/kernel:0dense_8/kernel/Assign$dense_8/kernel/Read/ReadVariableOp:0(2dense_8/random_uniform:08
^
dense_8/bias:0dense_8/bias/Assign"dense_8/bias/Read/ReadVariableOp:0(2dense_8/Const:08

Adam_1/iterations:0Adam_1/iterations/Assign'Adam_1/iterations/Read/ReadVariableOp:0(2-Adam_1/iterations/Initializer/initial_value:08

Adam_1/learning_rate:0Adam_1/learning_rate/Assign*Adam_1/learning_rate/Read/ReadVariableOp:0(20Adam_1/learning_rate/Initializer/initial_value:08
{
Adam_1/beta_1:0Adam_1/beta_1/Assign#Adam_1/beta_1/Read/ReadVariableOp:0(2)Adam_1/beta_1/Initializer/initial_value:08
{
Adam_1/beta_2:0Adam_1/beta_2/Assign#Adam_1/beta_2/Read/ReadVariableOp:0(2)Adam_1/beta_2/Initializer/initial_value:08
w
Adam_1/decay:0Adam_1/decay/Assign"Adam_1/decay/Read/ReadVariableOp:0(2(Adam_1/decay/Initializer/initial_value:08
I
	total_1:0total_1/Assigntotal_1/Read/ReadVariableOp:0(2	Const_3:08
I
	count_1:0count_1/Assigncount_1/Read/ReadVariableOp:0(2	Const_4:08
I
	total_2:0total_2/Assigntotal_2/Read/ReadVariableOp:0(2	Const_5:08
I
	count_2:0count_2/Assigncount_2/Read/ReadVariableOp:0(2	Const_6:08
I
	total_3:0total_3/Assigntotal_3/Read/ReadVariableOp:0(2	Const_7:08
I
	count_3:0count_3/Assigncount_3/Read/ReadVariableOp:0(2	Const_8:08
I
	total_4:0total_4/Assigntotal_4/Read/ReadVariableOp:0(2	Const_9:08
J
	count_4:0count_4/Assigncount_4/Read/ReadVariableOp:0(2
Const_10:08
J
	total_5:0total_5/Assigntotal_5/Read/ReadVariableOp:0(2
Const_11:08
J
	count_5:0count_5/Assigncount_5/Read/ReadVariableOp:0(2
Const_12:08* 
serving_default
2
default'
dense_1_input:0џџџџџџџџџ:
default/
lambda_1/l2_normalize:0џџџџџџџџџtensorflow/serving/predict*
applicability
4
default)
dropout_4_input:0џџџџџџџџџ9
default.
activation_8/Identity:0џџџџџџџџџtensorflow/serving/predict