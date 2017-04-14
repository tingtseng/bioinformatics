<%@ page contentType="text/html;charset=big5"%>
<%@ page import="java.sql.*"%>
<%@ page import="java.io.*"%>
<%@ page import="java.util.*"%>
<%@ page import="java.text.*"%>

<head>

	<link rel="icon" href="http://syslab.nchu.edu.tw/IIIDB/images/syslab.ico" type="image/x-icon" />
    <link rel="shortcut icon" href="http://syslab.nchu.edu.tw/IIIDB/images/syslab.ico" type="image/x-icon" />
	<link rel="bookmark" href="http://syslab.nchu.edu.tw/IIIDB/images/syslab.ico" type="image/x-icon" />

    <title>IIIDB - welcome</title>
    
    <meta http-equiv="content-type" content="application/xhtml+xml; charset=UTF-8" />
    
        <link rel="stylesheet" type="text/css" href="css/simple.css" media="screen, tv, projection" title="Default" />
	<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4/jquery.min.js"></script>
	<script>
		!window.jQuery && document.write('<script src="/jquery/jquery-1.4.3.min.js"><\/script>');
	</script>
	<script type="text/javascript" src="./jquery/fancybox/jquery.mousewheel-3.0.4.pack.js"></script>
	<script type="text/javascript" src="./jquery/fancybox/jquery.fancybox-1.3.4.pack.js"></script>
	<link rel="stylesheet" type="text/css" href="./jquery/fancybox/jquery.fancybox-1.3.4.css" media="screen" />

	<script type="text/javascript">
		$(document).ready(function() {
			$("a[rel=example_group]").fancybox({
				'transitionIn'		: 'none',
				'transitionOut'		: 'none',
				'titlePosition' 	: 'over',
				'titleFormat'		: function(title, currentArray, currentIndex, currentOpts) {
					return '<span id="fancybox-title-over">Image ' + (currentIndex + 1) + ' / ' + currentArray.length + (title.length ? ' &nbsp; ' + title : '') + '</span>';
				}
			});		
		});
	</script>
    

<script type="text/javascript" src="jquery/jquery.tablesorter.js"></script>
<script type="text/javascript">
	$(function() {
		$("table").tablesorter({widgets: ['zebra']});
	});
</script>
<style type="text/css">
body {
	background-image: url(./images/bg1.jpg);
	background-repeat: no-repeat;
}
</style>
</head>
<body id="your-site-id">

<div id="page" class="classname">
    
    <div id="header">
        <h1><a href="index.html">IIIDB</a></h1>
      <strong>A database for Isoform-Isoform Interactions and isoform network modules.</strong></div>
    
    <div id="wrapper"> 
        
        <div id="content">
        
            <div id="path">
                You are here: <a href="index.html">Home</a> 
                &raquo; <a href="Interact.html">Interaction Search</a>
            </div>
    
            <div id="main">

               




<%
request.setCharacterEncoding("big5");
String ID = request.getParameter("TranGeneID");
ID = ID.trim();
String action = request.getParameter("action");
String int_table = "";
if (action.equals("high")){
	int_table = "isoform_int_3";
	%>
     <h1>Interaction Search (High-confidence prediction, score > 3)</h1>
                <br />
    <%
}else{
	int_table = "isoform_int_2";
		%>
     <h1>Interaction Search (Low-confidence prediction, score > 2)</h1>
                <br />
    <%
}
Class.forName("com.mysql.jdbc.Driver");
Connection conn = DriverManager.getConnection("jdbc:mysql://140.120.203.50:3306/ting?useUnicode=true&characterEncoding=big5","ting", "ting2010");
Statement stmt = conn.createStatement();
Statement stmt2 = conn.createStatement();
NumberFormat f = NumberFormat.getNumberInstance();
f.setMaximumFractionDigits(5);

int i = 0;
String genesym = "null";
String geneid = "null";
String get_mod_sql="select genesym1, geneid1 from " + int_table + " WHERE geneid1= '" + ID + "' OR genesym1= '" + ID + "';" ;
if (stmt.execute(get_mod_sql))   {
   ResultSet rs = stmt.getResultSet();
   while (rs.next())
    {
		i++;
		genesym = rs.getString("genesym1");
		geneid = rs.getString("geneid1");
	}
   rs.close();
}
if (genesym.equals("null") && geneid.equals("null"))
{
%>	
	<h3>There are 0 interactions in your request: <%=ID%></h3>
<%
}else{
%>
	<h3>There are <%=i%> 
	interactions in your request: 
	<a href="http://www.ncbi.nlm.nih.gov/gene?term=<%=ID%>" target=_blank><%=genesym%> (GeneID: <%=geneid%>)</a>
	</h3>
<%



String sql="select * from " + int_table + " WHERE geneid1= '" + ID + "' OR genesym1= '" + ID + "';" ;
if (stmt.execute(sql))   {
	ResultSet rs = stmt.getResultSet();		
%>


<br />

<table cellspacing="0" class="tablesorter">
<thead>
    <TR>
    <th nowrap><B>Gene Symbol (Transcript ID)</B><span id="space">&nbsp;</span></th>
    <th nowrap><B>Gene ID</B><span id="space">&nbsp;</span></th>
    <th nowrap><B>Protein ID</B><span id="space">&nbsp;</span></th>
    <th width="43" nowrap><B>Score</B><span id="space">&nbsp;</span></th>
    <th width="33" nowrap><B>Co-exp</B><span id="space">&nbsp;</span></th>
    <th width="51" nowrap><B>Domain</B><span id="space">&nbsp;</span></th>
    <th width="30" nowrap><B>GO</B><span id="space">&nbsp;</span></th>
    <th width="65" nowrap><B>Ortholog</B><span id="space">&nbsp;</span></th>
    <th width="55" nowrap><B>Co-citation</B><span id="space">&nbsp;</span></th>
    </TR>   
</thead>
<tbody>
    <%

    while (rs.next())
    {
        String geneid1 = rs.getString("geneid1");
		String geneid2 = rs.getString("geneid2");
		String pid1 = rs.getString("proteinid1");
		if (pid1.equals("null")){
			pid1 = "N/A";
		}
		String pid2 = rs.getString("proteinid2");
		if (pid2.equals("null")){
			pid2 = "N/A";
		}
		%>
        <TR>
        <TD nowrap="nowrap"><B><a href="Interact.jsp?TranGeneID=<%= rs.getString("genesym1") %>&action=<%= action %>"><%= rs.getString("genesym1") %></a> (<%= rs.getString("tid1") %>) - <a href="Interact.jsp?TranGeneID=<%= rs.getString("genesym2") %>&action=<%= action %>"><%= rs.getString("genesym2") %></a> (<%= rs.getString("tid2") %>)</B></TD>
        <TD nowrap><B><%= geneid1 %> - <%= geneid2 %></B></TD>
        <TD nowrap><B><%= pid1 %> - <%= pid2 %></B></TD>
        <TD nowrap><B><%= f.format(rs.getDouble("score")) %></B></TD> 
        <TD nowrap><B><%= f.format(rs.getDouble("exp")) %></B></TD>
        <TD nowrap><B><%= f.format(rs.getDouble("domain")) %></B></TD>
        <TD nowrap><B><%= f.format(rs.getDouble("GO")) %></B></TD>
        <TD nowrap><B><%= f.format(rs.getDouble("ortholog")) %></B></TD> 
		<TD align="center"><B><%
		 //String co_str = "";
		 String co_sql="SELECT TA.pubmedid as puid FROM gene2pubmed2 AS TA, gene2pubmed2 AS TB	WHERE TA.geneid='" + geneid1 + "' AND TB.geneid='" + geneid2 + "'	AND TA.pubmedid=TB.pubmedid;" ;
		 if (stmt2.execute(co_sql)){
			ResultSet co_rs = stmt2.getResultSet();	
			if (co_rs.isAfterLast()!=co_rs.isBeforeFirst()){
				while (co_rs.next())
    			{
					%><a href="http://www.ncbi.nlm.nih.gov/pubmed?term=<%=co_rs.getString("puid")%>" target=_blank><%=co_rs.getString("puid")%></a><%
					out.print(" ");
			    	//co_str += co_rs.getString("puid") + " ";
				}
			}else{
		 		out.print("-");
				//co_str = "-";
		 	}
			co_rs.close();
		 }
		 %></B></TD>
		</TR>
        <%
    }
   
	rs.close();
	
%>
</tbody>
</TABLE>
<%
	
}
   
stmt.close();  
conn.close();  
}
%>  

<br>

<br>
                
               
                
            </div>

        </div>
        <div id="left">
          <div id="nav">
            <h3>Site Map</h3>
            <div class="inner_box">
              <ul id="mb73g0ebul_table" class="mb73g0ebul_menulist" style="width: 180px; height: 150px;">
                <li class="spaced_li"><a href="index.html#index"><img id="mbi_mb73g0_1" src="menu_images/mb_home.gif" name="mbi_mb73g0_1" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Home" title="" /></a></li>
                <li class="spaced_li"><a href="Interact.html"><img id="mbi_mb73g0_2" src="menu_images/mb_interaction_search.gif" name="mbi_mb73g0_2" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Interaction Search" title="" /></a></li>
                <li class="spaced_li"><a href="Enrich.html"><img id="mbi_mb73g0_3" src="menu_images/mb_module_search.gif" name="mbi_mb73g0_3" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Module Search" title="" /></a></li>
                <li class="spaced_li"><a href="index.html#about"><img id="mbi_mb73g0_5" src="menu_images/mb_about_us.gif" name="mbi_mb73g0_5" width="160" height="31" style="vertical-align: bottom;" border="0" alt="About Us" title="" /></a></li>
                <li class="spaced_li"><a href="index.html#tutorial"><img id="mbi_mb73g0_4" src="menu_images/mb_tutorial.gif" name="mbi_mb73g0_4" width="160" height="31" style="vertical-align: bottom;" border="0" alt="Tutorial" title="" /></a></li>
              </ul>
              <p>&nbsp;</p>
              <p>&nbsp;</p>
            </div>
            <link rel="stylesheet" href="menu_images/mbcsmb73g0.css" type="text/css" />
            <script type="text/javascript" src="menu_images/mbjsmb73g0.js"></script>
            <!-- i love clean source code -->
          </div>
          <div class="left_box">
            <h3>News and Updates</h3>
            <div class="inner_box">
              <p> <b>July 9, 2011:</b><br />
                New releases and related tools <br />
                <br />
                <a href="index.html">Read more...</a> </p>
            </div>
          </div>
          <div class="left_box">
            <h3>Related Tools</h3>
            <div class="inner_box">
              <p><a href="http://jigsaw.w3.org/css-validator/check/referer"> <img style="border:0;width:88px;height:31px"
        src="http://jigsaw.w3.org/css-validator/images/vcss-blue"
        alt="Valid CSS!" /> </a> </p>
<p> <a href="http://validator.w3.org/check?uri=referer"><img
      src="http://www.w3.org/Icons/valid-xhtml10" alt="Valid XHTML 1.0 Strict" height="31" width="88" /></a> </p>
            </div>
          </div>
        </div>
    </div>
    
    <div id="footer">
        <p>
        Copyright &copy; 2011 <a href="http://syslab.nchu.edu.tw/" target=_blank>Syslab</a>.</p>
    </div>

</div>
<a style="display:scroll;position:fixed;bottom:0px;right:5px;" href="#" title="" onFocus="if(this.blur)this.blur()"><img alt='' border='0' align="right" onMouseOver="this.src='http://4.bp.blogspot.com/-AKP9Gl_ets0/TmzrRBJPg_I/AAAAAAAADpY/M3KVvsDNJbA/s1600/Top.png'" src="http://1.bp.blogspot.com/-_gUa3K1wDNo/TmzrHROhBVI/AAAAAAAADpU/vOKm_7zL8DQ/s1600/Top_medium.png" onMouseOut="this.src='http://1.bp.blogspot.com/-_gUa3K1wDNo/TmzrHROhBVI/AAAAAAAADpU/vOKm_7zL8DQ/s1600/Top_medium.png'" /></a>
<!--BACKtoTOP-START--><!--BACKtoTOP-STOP-->
</body>
</html>
