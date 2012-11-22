//TO DO  - clean this whole thing up form jquery forms

var schema = new Object(); 

schema.create_new_field = function(token, value, permission) {
	var next_div_number = Math.floor(Math.random()*10000000);
	
	token = typeof token === 'undefined' ? '' : 'value="' +token +'"';
	value = typeof value === 'undefined' ? '' : 'value="' +value +'"';
	permission = typeof permission === 'undefined' ? 1 : 3;
	
	var div_text = '<div id="field' + next_div_number + '" class="schema_fields-container" uuid="' + next_div_number + '"> \
						<div class="span-6" schema="token"><input name="token" id=field' + next_div_number + '-1" type="text" ' + token + ' maxlength="255" /></div> \
						<div class="span-10" schema="value"><input name="value" id=field' + next_div_number + '-2" type="text" ' + value + ' maxlength="255" /></div> \
						<div class="span-5 last" schema="permission"> \
							<select name="permission" id="permission' + next_div_number + '" type="menu"> \
								<option ' + (permission == 1 ? "selected " : '') + 'value="1">Only Me</option>  \
								<option ' + (permission == 2 ? "selected " : '') + 'value="2">Visible on Domain</option> \
								<option ' + (permission == 3 ? "selected " : '') + 'value="3">Visible Anywhere</option> \
							</select> \
						\
						</div> \
						<div class="form-delete-row" schema_delete="fields"></div> \
					</div>';
	return div_text;
}

schema.create_new_editor = function(name, editables) {
	var next_div_number = Math.floor(Math.random()*10000000);
	
	name = typeof name === 'undefined' ? '' : 'value="' +name +'"';
	editables = typeof editables === 'undefined' ? [1,0,0,0] : editables; //editables is an array of booleans
	
	var div_text = '<div id="editor' + next_div_number + '" class="schema_editor-container" uuid="' + next_div_number + '"> \
						<div class="span-6" schema="name"><input name="name" type="text" ' + name + ' maxlength="255" /></div> \
						<div class="span-3 last" schema="editable1-outer"><input schema="editable1" type="checkbox" requires="" dependencies="2,3,4" class="rectCheck" ' + (editables[0] == 1 ? "checked " : '') + '/></div> \
						<div class="span-3" schema="editable2-outer"><input schema="editable2" type="checkbox" requires="1" dependencies="3,4" class="rectCheck" ' + (editables[1] == 1 ? "checked " : '') + '/></div> \
						<div class="span-3 last" schema="editable3-outer"><input schema="editable3" type="checkbox" requires="1,2" dependencies="4" class="rectCheck" ' + (editables[2] == 1 ? "checked " : '') + '/></div> \
						<div class="span-3 last" schema="editable4-outer"><input schema="editable4" type="checkbox" requires="1,2,3" dependencies="" class="rectCheck" ' + (editables[3] == 1 ? "checked " : '') + '/></div> \
						<div class="form-delete-row" schema_delete="editors"></div> \
					</div>';
	return div_text;
}

schema.populate_fields = function() {
	//this is the json containing all the data already present
	var demo_json = '{"user":"Susy Doozie", "fields": [ {"token": "firstname", "value": "Suzy", "permission": "3"}, {"token": "email", "value": "sdoozie@yahoo.com", "permission": "3"} ] }';
	var parsed = JSON.parse(demo_json);
			
	for (i = parsed.fields.length-1; i >= 0; i -= 1) {
		var temp_div = schema.create_new_field(parsed.fields[i].token, parsed.fields[i].value, parsed.fields[i].permission);
		$("#fields_container").prepend(temp_div);
	}
}

schema.populate_editors = function() {
	//this is the json containing all the data already present
	var demo_json = '{"editors": [ {"name": "Bob Jones", "editables": [1,1,0,0]}, {"name": "Sophia Werth", "editables": [0,1,1,1]} ] }';
	var parsed = JSON.parse(demo_json);
			
	for (i = parsed.editors.length - 1; i >= 0; i -= 1) {
		var temp_div = schema.create_new_editor(parsed.editors[i].name, parsed.editors[i].editables);
		$("#editors_container").prepend(temp_div);
	}
}

//note, this is a json object, but inside are just arrays
schema.undo_delete_jsons = JSON.parse('{"fields": [], "editors": []}');
//basically the delete function... hides a div and tells the server to update, but undo just re-shows it (easier than re-adding)
schema.activate_closes = function(activateMe) {
	$(activateMe).click(function() {
		var undo_div = $(this).parent().attr("uuid");
		undo_type = $(this).attr("schema_delete");

		//return json if possible, otherwise just save the HTML (in weird circumstances)
		var undo_content = $(this).parent().html();
		switch (undo_type) {
			case "fields":
				var token = $(this).siblings("div[schema=token]").children().val();
				var value = $(this).siblings("div[schema=value]").children().val();
				var permission = $(this).siblings("div[schema=permission]").children().val();
				var temp_undo_json = '{"jquery_div": "' + undo_div + '", "token": "' + token + '", "value": "' + value + '", "permission": "' + permission + '"}';
				schema.undo_delete_jsons.fields.push(temp_undo_json);
				$("#fields-undo-button").fadeIn("fast");
				break;
			case "editors": 
				var name = $(this).parent().children("div[schema=name]").children().val();
				var temp_permissions = [];
				$.each([1,2,3,4], function(index, value)  {
					temp_permissions.push($(activateMe).siblings("div[schema=editable"+value+"-outer]").children("input[type=checkbox]").is(":checked") ? 1 : 0);
				});
				var temp_undo_json = '{"jquery_div": "' + undo_div + '", "name": "' + name + '", "editables": [' + temp_permissions + '] }';
				schema.undo_delete_jsons.editors.push(temp_undo_json);
				$("#editors-undo-button").fadeIn("fast");
				break;
		}
		//hide it
		$("[uuid="+undo_div+"]").slideUp("slow"); ////TO DO - CREATE A METHOD THAT DOES THIS - DIV BY UUID
		
		//send to server
		/*
		libsend.call("remove"+undo_type, undo_div);
		*/
	});
}

//clicking on editor checkvoxes also activates/inactivates requirement/dependencies, respectively
schema.activate_editor_check = function(activateMe) {
	$(activateMe).click(function(e) {
		//e.preventDefault();
		var box = $(activateMe);
		if (box.is(":checked")) {
			var requirements = box.attr("requires").split(",");
			$.each(requirements, function(index, value)  {
				box.parent().siblings("div[schema=editable"+value+"-outer]").children("input[type=checkbox]").attr("checked", true);
			});
		} else {
			var dependencies = box.attr("dependencies").split(",");
			$.each(dependencies, function(index, value)  {
				box.parent().siblings("div[schema=editable"+value+"-outer]").children("input[type=checkbox]").attr("checked", false);
			});
			//don't want completely blank, reinstitute 1 if clicked
			if (box.attr("schema") == "editable1") { box.attr("checked", true) }						
		}	
	});
}

//must pass "field" or "editor"
schema.undo_delete = function(group) {
	if (group == "fields") {
		var temp_json = JSON.parse(schema.undo_delete_jsons.fields.pop());
		console.log('fields array reads after undo: ' + schema.undo_delete_jsons.fields +' .... ' + schema.undo_delete_jsons.fields[0]);
		//$("#fields_container").html(schema.create_new_field(temp_json.token, temp_json.value, temp_json.permission));
		$("[uuid="+temp_json.jquery_div+"]").slideDown("slow");
		if (!(schema.undo_delete_jsons.fields[0])) {
			$("#fields-undo-button").fadeOut("fast");
		}
	} else if (group == "editors") {
		var temp_json = JSON.parse(schema.undo_delete_jsons.editors.pop());
		console.log('editors array reads after undo: ' + schema.undo_delete_jsons.editors);
		//$("#editors_container").prepend(schema.create_new_editor(temp_json.name, temp_json.editables));
		$("[uuid="+temp_json.jquery_div+"]").slideDown("slow");
		if (!(schema.undo_delete_jsons.editors[0]) ){
			$("#editors-undo-button").fadeOut("fast");
		}
	} else {
	}
	//TO DO send to server
}

$(document).ready(function() {
	$(".schema_fields-buttons .form-add-row").parent().click(function() {
		$("#fields_container").append(schema.create_new_field()).children().last().hide().slideDown("slow"); 
		schema.activate_closes($("#fields_container").children().last().children(".form-delete-row")); //MAY WANT TO DO THIS A BETTER WAY
	});
	
	$(".schema_editors-buttons .form-add-row").parent().click(function() {
		$("#editors_container").append(schema.create_new_editor()).children().last().hide().slideDown("slow");
		schema.activate_closes($("#editors_container").children().last().children(".form-delete-row"));
	});
	
	$("#fields-undo-button").click(function(e) {
		e.preventDefault();
		schema.undo_delete("fields");
	});
	
	$("#editors-undo-button").click(function(e) {
		e.preventDefault();
		schema.undo_delete("editors");
	});
	
	$("input[name=schema_post-changes]").click(function(e) {
		e.preventDefault();
		//send update to server
	});
	
	$("input[name=schema_discard-changes]").click(function(e) {
		e.preventDefault();
		//send update to server
	});
	
	schema.populate_fields();
	schema.populate_editors();
	$(".form-delete-row").each(function() { schema.activate_closes(this) });
	$(".schema_editor-container input[type=checkbox]").each(function() { schema.activate_editor_check(this) });
	
});