/* Gets a list of images from the server and adds them to the imageList
 * Also adds a unique id to each part, viewable when clicked
 */
$(document).ready(function() {
    var deviceCount = 0;
    var command = {"command": "read"};
    $.get("EugeneServlet", command, function(response) {
        var toAppend = '<table class="table table-bordered table-hover" id="partsList"><thead><tr><th>Name</th><th>Type</th></tr></thead><tbody>';
        $.each(response["result"], function() {
            if (this["Type"].toLowerCase() !== "composite") {
                toAppend = toAppend + '<tr><td>' + this["Name"] + '</td><td>' + this["Type"] + '</td></tr>';
            }
        });
        toAppend = toAppend + "</tbody></table>";
        $('#partsListArea').html(toAppend);
        $("#partsList").dataTable({
            "bPaginate": false,
            "bScrollCollapse": true
        });
    });


    command = {"command": "imageList"};
    //render parts icons
    $.get("EugeneServlet", command, function(response) {
        var i = 0;
        $('#iconArea').html("");
        $.each(response, function() {
            var type = this["fileName"].split("\.")[0];
            $("#iconArea").append('<div class="span5"><li class="draggable partIcon" title= "' + type.replace(/-/g, ' ') + '" id="' + type + '"><div class="thumbnail"><img class="img-rounded" style="width:40px;height:80px" src="images/sbol_visual_jpeg/' + this["fileName"] + '"></div></li></div>');
            $('#' + type).dblclick(function() {
                //When clicked gives id name
                alert("creating new " + $(this).attr('id'));
            });
            i = i + 1;
        });
        $('#iconArea .draggable').draggable({
            helper: "clone",
            connectToSortable: ".sortable, #trash",
            revert: "invalid"
        });
        $('#iconArea li').on("click", function() {
            $(".selected").removeClass("selected");
            $(this).parent().addClass("selected");
            refreshPartsList($(this).attr("title"));
        });
    });


    //load parts list


    //Functions
    var refreshPartsList = function(s) {
        //TODO render parts list
        $('#partsList').html(s);
    };


    //Event Handlers
    $('#startButton').click(function() {
        $('#newDeviceButton').trigger("click");
    });

    $('#newDeviceButton').click(function() {
        if (deviceCount === 0) {
            $('#spectaclesEditorArea button').remove();
            $('#trash').droppable({
                hoverClass: "droppable-hover",
                drop: function(event, ui) {
                    var element = ui.draggable.css('position', '');
                    if (element.hasClass('blank')) {
                        element.parent().remove();
                    }
                    if (!element.hasClass("partIcon")) {
                        $(this).append(element);
                        $(ui.draggable).fadeOut(500);
                    }
                }}
            );
        }
        $('#spectaclesEditorArea').append('<ul id="device' + deviceCount + '" class="device droppable sortable"><li class="blank" style="vertical-align:bottom;height:80px;width:80px;border:1px solid grey" title="Drag a part here to get started">New Design</li></ul>');
        $("#device" + deviceCount).sortable({
            revert: true
        });
        $("#device" + deviceCount).droppable({
            drop: function() {
                $(this).droppable("destroy");
                var firstItem = $(this).find("li.blank");
                firstItem.html('<strong>' + $(this).attr("id") + '</strong>');
                firstItem.addClass("notSortable");
                $(this).sortable({
                    revert: true,
                    items: "li:not(.notSortable)",
                    connectWith: "#trash",
                    stop: function(event, ui) {
                        var firstItem = $(ui.item).parent().find("li.blank");
                        $(ui.item).parent().prepend(firstItem);
                        ui.item.removeClass('partIcon');
                    }
                });
                firstItem.attr("title", "Select to edit device properties");
                $(this).prepend(firstItem);
                firstItem.attr("style", 'height:40px;width:80px;background:#0081c2;vertical-align:top;padding:5px');
                firstItem.addClass("rotatedText");
                firstItem.click(function() {
                    $('.selected').removeClass("selected");
                    $(this).addClass("selected");
                });
            }
        });
        deviceCount = deviceCount + 1;
    });

    $('#runButton').click(function() {
        var command = {};
        command["command"] = "run";
        command["devices"] = "";
        command["parts"] = [];
        var devices = "";
        $('#spectaclesEditorArea ul').each(function() {
            if ($(this).find("li").length > 1) {
                var device = "Device " + $(this).attr("id") + "(";
                var count = 0;
                $(this).find("li").each(function() {
                    if (count > 0) {
                        device = device + $(this).attr("title").split(' ').join("") + ",";
                    }
                    count = count + 1;
                });
                device = device.substring(0, device.length - 1);
                device = device + ")";
                devices = devices + device + "|";
            }
        });
        devices = devices.substring(0, devices.length - 1);
        command["devices"] = devices;
        alert(command["devices"]);
        window.location.replace("http://cidar.bu.edu/ravencad/ravencad.html");
//        $.get("EugeneServlet", command, function(response) {
//            alert(response);
//        });
    });

});

