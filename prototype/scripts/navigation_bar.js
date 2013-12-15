$(document).ready(function(){
	$("#seeTabs").click(function(e) {
		e.preventDefault();
		if (!navigation_bar.showingTabs) {
			libsend.call("getTabList", "");
			navigation_bar.showingTabs = true;
		} else {
			navigation_bar.emptyTabNav();
		}
	});
});

var navigation_bar = new Object();

navigation_bar.showingTabs = false;

navigation_bar.activateTabNav = function() {
	$("li.openTab:not(.disabled)").click(function(e) {
		var socket_id = $(this).attr("socket_id");
		//if the url is empty it should just focus on the window with this uuid
		navigation_bar.emptyTabNav();
		libpage.add("EDITOR", socket_id);

	});
	
	$("li.createNewTab:not(.disabled)").click(function() {
		//interact with server appropriately
		//get a new uuid for the window as windowName

		navigation_bar.emptyTabNav();
		var newWindow = libpage.add("HOMEPAGE", "");
	});
	
	//set binding so mousing outside of tabs will close it after a delay
	var tabHideTimer;
	$("#tab_nav").on("mouseleave", (function(e) {
		tabHideTimer = setTimeout("navigation_bar.emptyTabNav()", 400);
	}));
	
	$("#tab_nav").on("mouseenter", (function(e) {
		clearTimeout(tabHideTimer);
	}));
};

navigation_bar.emptyTabNav = function() {
	$("#tab_nav").empty();
	navigation_bar.showingTabs = false;
};

navigation_bar.showTabs = function(tabs_list) {
	$("#tab_nav").empty();
	var tab_html = "";
	$.each(tabs_list, function(i, obj) {
		if (tabs_list[i].socket_id == window.name) {
			//make class disabled - <li  class="disabled">Tab Disabled</li>
			tab_html += '<li class="openTab disabled" socket_id="' + tabs_list[i].socket_id + '">' + tabs_list[i].name + ' -- ' +tabs_list[i].mode + '</li>';
		} else {
			tab_html += '<li class="openTab" socket_id="' + tabs_list[i].socket_id + '">' + tabs_list[i].name + ' -- ' +tabs_list[i].mode + '</li>';
		}
	});
	
	$("#tab_nav").append('<span class="arrow white topArrow"></span><ul>' + tab_html + '<li class="createNewTab">+ New Tab</li></ul>');
	navigation_bar.activateTabNav();
};
