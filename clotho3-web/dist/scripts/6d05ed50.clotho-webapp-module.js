angular.module("clotho.webapp", ["clotho.foundation", "clotho.interface", "clotho.dna", "ngSanitize", "ngRoute"]), angular.module("clotho.webapp").controller("HomeCtrl", ["$scope", "Clotho", "hotkeys", function(a, b, c) {
  a.modalContent = "<p>Welcome to Clotho!</p><p>Clotho is a platform for automating your genetic engineering projects. Learn how to use Clotho by starting the trail below!</p>", a.downloadClotho = function() {
    window.open("https://github.com/CIDARLAB/clotho3", "_blank")
  }, a.browseApps = function() {
    window.open('/#!/apps');
    //window.open("https://github.com/CIDARLAB/clotho3apps", "_blank")
  }, a.enterClotho = function() {
    window.open("http://www.synbiotrails.org", "_blank")

    //b.startTrail("org.clothocad.trails.LearningClotho")
  }, c.bindTo(a).add({
    combo: "h",
    description: "Show Intro Modal",
    callback: function() {
      a.showHelp = !a.showHelp
    }
  })
}]), angular.module("clotho.webapp").controller("SettingsCtrl", ["$scope", "Clotho", function(a, b) {
  b.get("clotho.developer.maxbates", {
    mute: !0
  }).then(function(b) {
    a.person = b
  })
}]), angular.module("clotho.webapp").controller("TeamCtrl", ["$scope", "Clotho", function(a, b) {
  b.query({
    id: {
      $regex: "clotho.developer.*",
      $options: "i"
    }
  }, {
    mute: !0
  }).then(function(b) {
    a.ClothoTeam = b
  })
}]), angular.module("clotho.webapp").controller("BrowserCtrl", ["$scope", "Clotho", "$filter", "ClothoSchemas", function(a, b, c, d) {
  a.orderers = [{
    name: "Name",
    criteria: "name",
    "class": "glyphicon-sort-by-alphabet"
  }, {
    name: "type",
    criteria: d.dirtyDetermineType,
    "class": "glyphicon-th-list"
  }, {
    name: "Schema",
    criteria: "schema",
    "class": "glyphicon-cog"
  }], a.setOrder = function(b) {
    a.currentOrder = b === a.currentOrder ? null : b
  }, a.filters = [{
    name: "Has Description",
    filter: function(a) {
      return angular.isDefined(a.description) && a.description
    },
    "class": "glyphicon glyphicon-comment"
  }], a.currentFilter = function() {
    return !0
  }, a.setFilter = function(b) {
    a.currentFilter = b === a.currentFilter ? function() {
      return !0
    } : b
  }, a.collections = [{
    name: "My Collection",
    author: "uniqueUserID",
    description: "my Collection of all the things.",
    items: {
      "clotho.developer.maxbates": {
        note1: "blah blah blah blah"
      },
      "clotho.enzyme.BglII": {
        note1: "yad yad ayayayayy"
      },
      "clotho.part.jtk2134": {
        note2: "bling bling blang"
      }
    }
  }], a.setCollection = function(c) {
    a.currentQuery != c && (a.currentQuery = c, a.resultArray = [], angular.forEach(c.items, function(c, d) {
      b.get(d, {
        mute: !0
      }).then(function(b) {
        a.resultArray.push(b)
      })
    }))
  }, a.queries = [{
    name: "Schemas",
    query: {
      schema: "org.clothocad.core.schema.Schema"
    }
  }, {
    name: "Parts",
    query: {
      schema: "org.clothocad.model.Part"
    }
  }, {
    name: "Vectors",
    query: {
      schema: "Vector"
    }
  }, {
    name: "People",
    query: {
      schema: "org.clothocad.model.LabPerson"
    }
  }], a.newQuery = {}, a.saveNewQuery = function() {
    a.queries.push(a.newQuery), a.newQuery = {}
  }, a.setCurrentQuery = function(c, d) {
    a.currentQuery != c && (d = d || 200, a.currentQuery = c, b.query(c, {
      maxResults: d
    }).then(function(b) {
      a.resultArray = b
    }))
  }, a.sort = function(b) {
    if (b) {
      if (a.catSort) return;
      a.catSort = !0, a.recent = c("categorize")(a.recent_array, "type"), a.recent.Instance = c("categorize")(a.recent.Instance, "schema.name")
    } else {
      if (!a.catSort && angular.isDefined(a.catSort)) return;
      a.catSort = !1, a.recent = {
        entries: angular.copy(a.recent_array)
      }
    }
  }, a.collectionIconClass = "glyphicon glyphicon-briefcase", a.filterIconClass = "glyphicon glyphicon-filter", a.queryIconClass = "glyphicon glyphicon-wrench", a.setCurrentQuery(a.queries[0].query)
}]), angular.module("clotho.webapp").controller("EditorCtrl", ["$scope", "$route", "$location", "Clotho", "ClothoSchemas", function(a, b, c, d, e) {
  a.$watch("editable.id", function(a) {
    c.search("id", a || null).replace()
  }), a.$on("$routeUpdate", function(b, c) {
    var d = c.params.id;
    angular.isEmpty(a.editable) || a.editable.id == d || (a.editableId = d)
  }), a.$on("$destroy", function() {
    c.search("id", null).replace()
  }), b.current.params.id && (a.editableId = b.current.params.id), a.editModePass = !1, a.objectTypes = e.sharableTypes, a.schemas = [], e.retrievedSchemas.then(function(b) {
    a.schemas = b
  }), a.queryWrapper = function(a) {
    return d.autocomplete(a).then(function(a) {
      return a || []
    })
  }, a.createNewNonInstance = function(b) {
    a.editable = e.createScaffold(b), a.editModePass = !0
  }, a.createNewInstance = function(b) {
    a.editable = {
      id: b.id
    }, a.editModePass = !0
  }, a.editExisting = function(b) {
    a.editableId = b.id, a.editModePass = !0
  }
}]), angular.module("clotho.webapp").controller("ExecutorCtrl", ["$scope", "$route", "$location", "$filter", "Clotho", "ClothoSchemas", function(a, b, c, d, e, f) {
  a.$on("$routeUpdate", function(b, c) {
    var d = c.params.id;
    a.functionId != d && (a.functionId = d)
  }), a.$on("$destroy", function() {
    c.search("id", null).replace()
  }), b.current.params.id && (a.functionId = b.current.params.id), a.$watch("functionId", function(b) {
    e.get(b).then(function(b) {
      f.isFunction(b) && (a.function = b)
    })
  }), a.$on("$routeUpdate", function(b, c) {
    var d = c.params.id;
    angular.isEmpty(a.editable) || a.editable.id == d || (a.functionId = d)
  }), a.$watch("function.id", function(a) {
    c.search("id", a || null).replace()
  }), a.queryWrapper = function(a, b) {
    return e.autocomplete(a).then(function(a) {
      return angular.isUndefined(b) ? a : d("filter")(a, function(a) {
        return "function" == b ? f.isFunction(a) : f.isInstanceOfSchema(a, b)
      })
    })
  }, a.onExecute = function() {}, a.editExisting = function(b) {
    f.isFunction(b) && (a.function = b)
  }
}]), angular.module("clotho.trails").controller("TrailCtrl", ["$scope", "$route", "$timeout", "Clotho", "Trails", "hotkeys", "$location", function(a, b, c, d, e, f, g) {
  a.id = b.current.params.id, a.trail = b.current.locals.trail, b.current.params.position && (a.current = b.current.params.position, a.currentPage = e.extractPage(a.trail, a.current)), a.activate = e.activate, a.favorite = function() {
    e.favorite(a.id)
  }, a.share = function() {
    e.share(a.id)
  }, a.next = function() {
    console.log("activating next", a.current, e.calcNextPage(a.trail, a.current)), e.activate(e.calcNextPage(a.trail, a.current))
  }, a.prev = function() {
    e.activate(e.calcPrevPage(a.trail, a.current))
  }, a.mapIcon = e.mapIcon, a.$on("$routeUpdate", function(b, c) {
    a.current != c.params.position && (a.current = c.params.position || null, a.currentPage = e.extractPage(a.trail, a.current), angular.isEmpty(a.currentPage) && a.activate("0-0"))
  }), a.$on("$destroy", function() {
    g.search("position", null).replace(), g.search("id", null).replace()
  }), a.activate(a.current)
}]), angular.module("clotho.webapp").controller("TrailsCtrl", ["$scope", "$location", "Clotho", function(a, b, c) {
  a.headerMockTrail = {
    name: "Trail Browser",
    description: "Learn Synthetic Biology with Clotho"
  }, a.trails = [], c.query({
    schema: "org.clothocad.model.Trail"
  }, {
    mute: !0
  }).then(function(b) {
    a.trails = b
  }), a.defaultTrailIcon = "images/trails/trails_logo.png", a.startTrail = function(a) {
    c.startTrail(a)
  }, a.startTrailPage = function(a, d) {
    b.search("position", d), c.startTrail(a)
  }, a.highlight = function(b) {
    a.highlighted = b, a.loading = !0, c.get(b.id, {
      mute: !0
    }).then(function(b) {
      a.loading = !1, a.selected = b
    })
  }
}]), angular.module("clotho.trails").controller("TrailSplashCtrl", ["$scope", "$location", "$http", "Clotho", function(a, b, c, d) {
  c.get("models/trail-splash.json", {
    cache: !0
  }).success(function(b) {
    a.topics = b
  }), a.startTrail = function(a) {
    d.startTrail(a)
  }, a.startTrailPage = function(a, c) {
    b.search("position", c), d.startTrail(a)
  }, a.highlight = function(b) {
    a.highlighted = b, a.loading = !0, d.get(b.id, {
      mute: !0
    }).then(function(b) {
      a.loading = !1, a.selected = b
    })
  }
}]), angular.module("clotho.webapp").controller("WidgetsCtrl", ["$scope", "ClientAPI", function(a, b) {
  a.myObj = {
    myProp: "aacccggttt"
  }, a.myModel = "WHATS up", a.myDNAmodel = "aaaagt", a.bootstrapCallback = function(a) {
    console.log("widget controller callback! passed element:", a)
  }, a.bootstrapNewApp = function() {
    b.display("123456789", "insert-widgets-here")
  }, a.bootstrapNewAppExternal = function() {
    b.display("123456789")
  }, a.bootstrapSimple = function() {
    b.display("813579135", "insert-widgets-here")
  }, a.someValue = "from parent controller not passed down", a.openCallback = function() {
    console.log("modal openeed")
  }, a.closeCallback = function() {
    console.log("modal closed")
  }
}]), angular.module("clotho.webapp").controller("ImportCtrl", ["$scope", "$location", "Clotho", function(a, b) {
  a.importable = [{
    name: "Youtube Playlist",
    route: "/import/youtubePlaylist"
  }, {
    name: "ApE Files",
    route: "/import/ape"
  }, {
    name: "NCBI Entrez Gene",
    route: "/import/ncbi"
  }], a.goToImportable = function(a) {
    b.path(a.route)
  }
}]), angular.module("clotho.webapp").controller("YoutubePlaylistImportCtrl", ["$scope", "Youtube", "Clotho", "$q", function(a, b, c) {
  a.playlistId = "PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy", a.$watch("playlistId", function(c) {
    b.playlistItems(c).then(function(b) {
      a.playlistInfo = b
    }), b.playlistInfo(c).then(function(b) {
      a.playlistItems = b
    }), b.playlistToTrail(c).then(function(b) {
      a.playlistTrail = b
    })
  }), a.create = function() {
    angular.isEmpty(a.playlistTrail) || c.create(a.playlistTrail)
  }, a.$watch("search", function(c) {
    c ? b.videoSearch(c).then(function(b) {
      a.searchResult = b
    }) : a.searchResult = ""
  })
}]), angular.module("clotho.webapp").controller("ApeImportCtrl", ["$scope", "Clotho", function(a, b) {
  a.converted = [], a.processGenbankFile = function(c, d) {
    b.run("org.andersonlab.py_convertGB", [d]).then(function(b) {
      console.log(b), a.converted.push(b[0])
    })
  }, a.createNucSeq = b.create
}]), angular.module("clotho.webapp").service("ImportCSVService", ["$window", "$q", "ClothoExtensions", function(a, b, c) {
  var d = this;
  d.ready = c.mixin("lib/babyparse.js").then(function() {
    return console.log(a.Baby), a.Baby
  }), d.defaults = {
    header: !0
  }, d.parse = function(a, c) {
    return angular.isString(a) ? (c = angular.extend({}, d.defaults, c), d.ready.then(function(b) {
      return b.parse(a, c)
    })) : b.when({})
  }
}]).controller("ImportCSVCtrl", ["$scope", "$q", "Clotho", "ImportCSVService", function(a, b, c, d) {
  function e() {
    d.parse.apply(null, arguments).then(function(b) {
      a.parsed = b, a.parsedData = b.data, b.data.length && (a.fieldDemo = b.data[0], a.fieldMap = {}, angular.forEach(Object.keys(a.fieldDemo), function(b) {
        a.fieldMap[b] = b
      }))
    })
  }
  a.csvText = "", a.options = {
    header: !0
  }, a.process = function() {
    angular.isDefined(a.csvText) && a.csvText.length && e(a.csvText, a.options)
  }, a.removeField = function(b) {
    delete a.fieldMap[b]
  }, a.reparse = function() {
    var b = [];
    angular.forEach(a.parsedData, function(c) {
      var d = {};
      angular.forEach(a.fieldMap, function(a, b) {
        d[a] = c[b]
      }), b.push(d)
    }), a.parsedData = b
  }, a.importParsed = function() {
    a.createdIds = [], angular.forEach(a.parsedData, function(b) {
      c.create(b).then(function(b) {
        a.createdIds.push(b)
      })
    })
  }
}]);
