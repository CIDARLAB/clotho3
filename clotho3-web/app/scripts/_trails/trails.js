angular.module('clotho.trails').service('Trails', function(Clotho, $q) {

	var compile = function TrailCompile(trail) {

		//If pass by reference (depending on Clotho.get() ) need to copy to don't edit in dependencies
		//trail = angular.copy(trail);

		var transcludes = trail.dependencies || null,
			deferred = $q.defer();

		if (!transcludes && !trail.mixin) {
			deferred.resolve(trail);
			return deferred.promise;
		}

		var final_contents = [],
			promises = [];

		//get the transcluded trails... will be fast if in collector already
		(transcludes) && angular.forEach(transcludes, function(id) {
			promises.push(Clotho.get(id));
		});

		$clotho.extensions.mixin(trail.mixin)
			.then(function() {
				return $q.all(promises)
			})
			//after download all, pluck out the chapters we need
			.then(function (downloads) {

				//reorganize transcludes so can reference by id
				transcludes = {};
				angular.forEach(downloads, function(transclude) {
					transcludes[transclude.id] = transclude;
				});

				//iterate through trail, pushing in chapters
				angular.forEach(trail.contents, function (chapter, ind) {

					if (typeof chapter.transclude == 'undefined') {
						final_contents.push(chapter);
					} else {
						//chapters to include :
						var chapterID = chapter.transclude.id,
							chapterNum = chapter.transclude.chapters;

						if ((chapterNum == "all") || (typeof chapterNum == 'undefined')) {
							for (var i = 0; i < transcludes[chapterID]['contents'].length; i++) {
								final_contents.push(transcludes[chapterID]['contents'][i]);
							}
						} else {
							var startStop = chapterNum.split("-");
							if (startStop.length == 1) {
								final_contents.push(transcludes[chapterID]['contents'][startStop[0]]);
							} else {
								if (startStop[0] > startStop[1])
									return "wrong format - start must be smaller than end";

								for (var i = startStop[0]; i <= startStop[1]; i++) {
									final_contents.push(transcludes[chapterID]['contents'][i])
								}
							}
						}
					}
				});

				trail.contents = final_contents;
				deferred.resolve(trail);
			});

		return deferred.promise;
	};

	//in form <Chapter>-<Page>
	var extractPage = function(Trail, indices) {
		var pos = indices.split("-");
		var page = Trail.contents[pos[0]]['pages'][pos[1]];
		return page;
	};

	//bring back in logic from trail-module.js
	var calcNextPage = function(Trail, oldpos) {
		oldpos = (typeof oldpos != 'undefined') ? oldpos.split("-") : [0, -1];
		var newpos;

		if (typeof Trail.contents[oldpos[0]]['pages'][+oldpos[1] + 1] != 'undefined')
			newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
		else if (typeof Trail.contents[+oldpos[0] + 1]['pages'] != 'undefined')
			newpos = (+oldpos[0] + 1) + '-' + 0;
		else {
			return;
		}
		return newpos;
	};

	var calcPrevPage = function(Trail, oldpos) {
		if (oldpos == '0-0') return;

		oldpos = (typeof oldpos != 'undefined') ? oldpos.split("-") : [0, 1];
		var newpos;

		if (typeof Trail.contents[oldpos[0]]['pages'][+oldpos[1] - 1] != 'undefined')
			newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
		else if (typeof Trail.contents[+oldpos[0] - 1]['pages'])
			newpos = (+oldpos[0] - 1) + '-' + (Trail.contents[+oldpos[0] - 1]['pages'].length - 1);
		else {
			return;
		}
		return newpos;
	};

	var trailIconMap = {
		'quiz'        : 'icon-pencil'     ,
		'list'        : 'icon-list-alt'   ,
		'eye'         : 'icon-eye-open'   ,
		'info'        : 'icon-info-sign'  ,
		'video'       : 'icon-film'       ,
		'template'    : 'icon-book'       ,
		'exercise'    : 'icon-edit'       ,
		'undefined'   : 'icon-file'         //fallthrough
	};

	var mapIcon = function(iconName) {
		iconName = iconName || 'undefined';
		return trailIconMap[iconName];
	};



	var favorite = function(id) {
		console.log("favorite trail with id: " + id);
	};


	return {
		compile : compile,
		extractPage : extractPage,
		calcNextPage : calcNextPage,
		calcPrevPage : calcPrevPage,
		mapIcon : mapIcon,
		share : Clotho.share,
		favorite : favorite
	}
});