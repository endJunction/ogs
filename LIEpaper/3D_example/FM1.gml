<?xml version="1.0"encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>TaskB_3D</name>
	<points>
		<point id="0"x="10"y="-10"z="-10"/>
		<point id="1"x="10"y="10"z="-10"/>
		<point id="2"x="10"y="10"z="10"/>
		<point id="3"x="10"y="-10"z="10"/>
		<point id="4"x="10"y="-4.6630766" z="-10"/>
		<point id="5"x="10"y="4.6630766" z="10"/>
		<point id="6"x="0"y="0"z="0"name="INJECTION_P0"/>
		<point id="7"x="0"y="-10"z="-10"/>
		<point id="8"x="0"y="10"z="-10"/>
		<point id="9"x="0"y="10"z="10"/>
		<point id="10"x="0"y="-10"z="10"/>
		<point id="11"x="0"y="-4.6630766" z="-10"/>
		<point id="12"x="0"y="4.6630766" z="10"/>
		<point id="13"x="0"y="-10"z="0"/>
		<point id="14"x="10"y="-10"z="0"/>
		<point id="15"x="0"y="10"z="0"/>
		<point id="16"x="10"y="10"z="0"/>
		<point	id="17"	x="10"	y="0"	z="0"/>
		<point	id="18"	x="0"	y="0"	z="0.25"	name="anchor_Pup"	/>
		<point	id="19"	x="0"	y="0"	z="-0.25"	name="anchor_Pdown"	/>
		<point	id="20"	x="0.0"	y="-0.633927393"z="-1.359461681"	name="P3"/>
		<point	id="21"	x="0.0"	y="0.633927393"z="1.359461681"	name="P3_MIRROR"/>
		<point	id="22"	x="1.5"	y="0.0"z="0.0"	name="M1"/>
		<point	id="23"	x="0.0"	y="0.023234652"z="0.049826873"	name="INJECTION_P5"/>
		<point	id="24"	x="0.0"	y="-0.023234652"z="-0.049826873"	name="INJECTION_P1"/>
		<point	id="25"	x="0.055"	y=".023234652"z=".049826873"	name="INJECTION_P4"/>
		<point	id="26"	x="0.055"	y="0.0"z="0.0"	name="INJECTION_P3"/>
		<point	id="27"	x="0.055"	y="-0.023234652"z="-.049826873"	name="INJECTION_P2"/>
	</points>
	<polylines>
		<polyline id="0"name="PLY_TOP">
			<pnt>10</pnt>
			<pnt>12</pnt>
			<pnt>9</pnt>
			<pnt>2</pnt>
			<pnt>5</pnt>
			<pnt>3</pnt>
			<pnt>10</pnt>
		</polyline>
		<polyline id="1"name="PLY_BOTTOM">
			<pnt>7</pnt>
			<pnt>11</pnt>
			<pnt>8</pnt>
			<pnt>1</pnt>
			<pnt>4</pnt>
			<pnt>0</pnt>
			<pnt>7</pnt>
		</polyline>
		<polyline id="2"name="PLY_LEFT">
			<pnt>3</pnt>
			<pnt>10</pnt>
			<pnt>7</pnt>
			<pnt>0</pnt>
			<pnt>3</pnt>
		</polyline>
		<polyline id="3"name="PLY_RIGHT">
			<pnt>2</pnt>
			<pnt>9</pnt>
			<pnt>8</pnt>
			<pnt>1</pnt>
			<pnt>2</pnt>
		</polyline>
		<polyline id="4"name="PLY_FRONT">
			<pnt>3</pnt>
			<pnt>5</pnt>
			<pnt>2</pnt>
			<pnt>1</pnt>
			<pnt>4</pnt>
			<pnt>0</pnt>
			<pnt>3</pnt>
		</polyline>
		<polyline id="5"name="PLY_BACK">
			<pnt>10</pnt>
			<pnt>12</pnt>
			<pnt>9</pnt>
			<pnt>8</pnt>
			<pnt>11</pnt>
			<pnt>7</pnt>
			<pnt>10</pnt>
		</polyline>
		<polyline id="6"name="PLY_IN">
			<pnt>13</pnt>
			<pnt>14</pnt>
		</polyline>
		<polyline id="7"name="PLY_OUT">
			<pnt>15</pnt>
			<pnt>16</pnt>
		</polyline>
		<polyline id="8"name="PLY_INJECTION">
			<pnt>6</pnt>
			<pnt>24</pnt>
			<pnt>27</pnt>
			<pnt>26</pnt>
			<pnt>25</pnt>
			<pnt>23</pnt>
			<pnt>6</pnt>
		</polyline>
		<polyline id="9"name="PLY_FTOP">
			<pnt>5</pnt>
			<pnt>12</pnt>
		</polyline>		
		<polyline id="10"name="PLY_FFRONT">
			<pnt>4</pnt>
			<pnt>5</pnt>
		</polyline>				
		<polyline id="11"name="PLY_FBOTTOM">
			<pnt>11</pnt>
			<pnt>4</pnt>
		</polyline>		
		<polyline id="12"name="PLY_FBACK">
			<pnt>11</pnt>
			<pnt>12</pnt>
		</polyline>				
	</polylines>
	<surfaces>
        <surface id="0"name="left"><!-- y=0 -->
            <element p1="3"p2="10"p3="7"/>
            <element p1="3"p2="7"p3="0"/>
        </surface>
        <surface id="1"name="right"><!-- y=1 -->
            <element p1="2"p2="9"p3="8"/>
            <element p1="2"p2="8"p3="1"/>
        </surface>
        <surface id="2"name="top"><!-- z=1 -->
            <element p1="10"p2="9"p3="2"/>
            <element p1="10"p2="2"p3="3"/>
        </surface>
        <surface id="3"name="bottom"><!-- z=0 -->
            <element p1="7"p2="8"p3="1"/>
            <element p1="7"p2="1"p3="0"/>
        </surface>
        <surface id="4"name="front"><!-- x=10 -->
            <element p1="3"p2="2"p3="1"/>
            <element p1="3"p2="1"p3="0"/>
        </surface>
        <surface id="5"name="back"><!-- x=0 -->
            <element p1="10"p2="9"p3="8"/>
            <element p1="10"p2="8"p3="7"/>
        </surface>
        <surface id="6"name="injection"> <!-- INJECTION AREA -->
            <element p1="23"p2="24"p3="27"/>
            <element p1="23"p2="27"p3="25"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
