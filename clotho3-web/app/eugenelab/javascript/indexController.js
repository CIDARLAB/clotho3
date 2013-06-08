//I want my javascript to run only once the page is rendered, this helps the page load faster
$(document).ready(function() {
    
    //javascript to control the drag and drop elements
    $(function() {
        $("#draggable").draggable({snap: true});
        $("#draggable2").draggable({snap: ".ui-widget-header"});
        $("#draggable3").draggable({snap: ".ui-widget-header", snapMode: "outer"});
        $("#draggable4").draggable({grid: [20, 20]});
        $("#draggable5").draggable({grid: [80, 80]});
    });
    //this is an event handler for the element with the id bigButton.
//    this handler decides what happens when bigButton is clicked
    $('#bigButton').click(function(){
        var data = {"command": "test"}; //this is the JSON object that I'm sending to DemoServlet
        //I'm going to do a get request from DemoServlet
        $.get("DemoServlet", data , function(response) { //i'm defining a callback function for the get request to the servlet DemoServlet
            alert(response); //I'm going to display the response of the servlet in an alert window
            //I'm also going to insert the response into my page
            $('#sandbox').html("<h1>"+response+"</h1>"); //im inserting the response into the element with id sandbox
        });
    });
});


