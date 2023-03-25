/*
 * Copyright 2022 University of California, Riverside
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.ucr.cs.bdlab.davinci

import edu.ucr.cs.bdlab.beast.geolite.{Feature, GeometryHelper, GeometryReader, IFeature}
import edu.ucr.cs.bdlab.beast.sql.CreationFunctions.geometryFactory
import edu.ucr.cs.bdlab.beast.util.BitArray
import edu.ucr.cs.bdlab.davinci.IntermediateTile.maxPointsPerTile
import org.geotools.geometry.jts.JTSFactoryFinder
import org.locationtech.jts.geom.{Coordinate, CoordinateSequence, CoordinateXY, Envelope, GeometryFactory, LinearRing, MultiPolygon, Polygon}

import scala.collection.mutable.ArrayBuffer

/**
 * A class that stores intermediate data to build vector tiles.
 *
 * @param resolution the resolution of the tile in pixels.
 */
class IntermediateTile(resolution: Int) extends Serializable {

  /** A list of features that are stored in full */
  val features: collection.mutable.ArrayBuffer[IFeature] = new collection.mutable.ArrayBuffer[IFeature]()

  /** Counts total number of points of all features in the list */
  var numPoints: Int = 0

  /** A bit map that marks which pixels are occupied with rasterized features */
  var occupiedPixels: BitArray = _

  var bufferSize = 0

  // TODO store additional aggregated attribute

  /**
   * Add the given feature to this tile.
   * The feature is either added as-is, or rasterized and aggregated, depending on the state of this tile.
   *
   * @param feature the feature to add
   * @return this tile to allow chaining other operations
   */
  def addFeature(feature: IFeature): IntermediateTile = {
    if (!this.isRasterized) {
      // Still collecting features in full
      features.append(feature)
      numPoints += feature.getGeometry.getNumPoints

      if (numPoints > maxPointsPerTile) {
        occupiedPixels = new BitArray(resolution * resolution)
        for (f <- features)
          rasterizeFeature(f)
        features.clear()
      }
    } else {
      rasterizeFeature(feature)
    }
    this
  }

  /** Whether this tile is in the rasterized space or not */
  def isRasterized: Boolean = this.occupiedPixels != null

  private def getPixelOffset(x: Int, y: Int): Long = {
    if (x < -bufferSize || x >= resolution + bufferSize || y < -bufferSize || y >= resolution + bufferSize) return -1
    return (y + bufferSize) * (resolution + 2 * bufferSize) + x
  }
  private def rasterizeFeature(feature: IFeature): Unit = {
    // rasterize the feature into occupied pixels
    val c: Coordinate=feature.getGeometry.getCoordinate
    val x: Int = (c.x.toInt)
    val y: Int = (c.y.toInt)
    val pixelOffset: Long = getPixelOffset(x, y)
    if (pixelOffset != -(1)) {
      occupiedPixels.set(pixelOffset, true)
    }

    // TODO aggregate additional attributes when rasterizing features
    val additionalAttr = new ArrayBuffer[Polygon]()
  }

  /**
   * Merge another tile into this tile
   *
   * @param other the other tile to merge into this one
   * @return this tile after merging with the other one
   */
  def merge(other: IntermediateTile): IntermediateTile = {
    if (this.isRasterized && other.isRasterized) {
      // Both are rasterized, just combine the two sets of occupied pixels
      this.occupiedPixels.inplaceOr(other.occupiedPixels)
    } else if (this.isRasterized && !other.isRasterized) {
      // Rasterize all incoming features into this
      for (f <- other.features)
        rasterizeFeature(f)
    } else if (!this.isRasterized && other.isRasterized) {
      // Rasterize current one and then merge with the other one
      occupiedPixels = new BitArray(resolution * resolution)
      for (f <- features)
        rasterizeFeature(f)
      features.clear()
      this.occupiedPixels.inplaceOr(other.occupiedPixels)
    } else {
      // None is rasterized. Add incoming features one-by-one and rasterize if needed
      for (f <- other.features)
        addFeature(f)
    }
    // Return this tile to chain merge operations if needed
    this
  }



  /**
   * Converts this intermediate tile into a final [[VectorTile.Tile]]
   * @return the vector tile that represents all features in this tile
   */
  def vectorTile: VectorTile.Tile = {
    val layer = if (!this.isRasterized) {
      val vectorLayerBuilder = new VectorLayerBuilder(resolution, "features")
      // Add all features one-by-one
      for (f <- features)
        vectorLayerBuilder.addFeature(f)
      vectorLayerBuilder.build()
    } else {
      // Add occupied pixels
      val vectorLayerBuilder = new VectorLayerBuilder(resolution, "aggregate")
      // The object delineation algorithm to aggregate occupied pixels
      val corners = new ArrayBuffer[Node]()
      var leftVertex: Node = null
      val topVertices: Array[Node] = new Array[Node](resolution + 1)
      var offset = 0

      for (y <- 0 to resolution) {
        for (x <- 0 to resolution) {
          val blocked0: Int = if (x > -bufferSize && y > -bufferSize && occupiedPixels.get(offset - resolution - 1)) 1 else 0
          val blocked1: Int = if (x < resolution + bufferSize && y > -bufferSize && occupiedPixels.get(offset - resolution)) 2 else 0
          val blocked2: Int = if (x > -bufferSize && y < resolution + bufferSize && occupiedPixels.get(offset - 1)) 4 else 0
          val blocked3: Int = if (x < resolution + bufferSize && y < resolution + bufferSize && occupiedPixels.get(offset)) 8 else 0
          val pixelType = blocked0 + blocked1 + blocked2 + blocked3
          pixelType match {
            case 0 | 3 | 5 | 10 | 12 | 15 => // Do nothing
            case 1 =>
              val newVertex = Node(x, y, topVertices(x))
              leftVertex.next = newVertex

            case 2 =>
              val newVertex = Node(x, y, null)
              topVertices(x).next = newVertex
              leftVertex = newVertex

            case 4 =>
              val newVertex = Node(x, y, leftVertex)
              topVertices(x) = newVertex

            case 6 =>
              val newVertex1 = Node(x, y, leftVertex)
              val newVertex2 = Node(x, y, null)
              topVertices(x).next = newVertex2
              topVertices(x) = newVertex2
              leftVertex = newVertex2
              topVertices(x) = newVertex1
            case 7 =>
              val newVertex = Node(x, y, null)
              topVertices(x) = newVertex
              leftVertex = newVertex
              corners.append(newVertex)
            case 8 =>
              val newVertex = Node(x, y, null)
              topVertices(x) = newVertex
              leftVertex = newVertex
              corners.append(newVertex)

            case 9 =>
              val newVertex1 = Node(x, y, topVertices(x))
              leftVertex.next = newVertex1
              val newVertex2 = Node(x, y, null)
              leftVertex = newVertex2
              topVertices(x) = newVertex2
              corners.append(newVertex2)

            case 11 =>
              val newVertex = Node(x, y, null)
              leftVertex.next = newVertex
              leftVertex = null
              topVertices(x) = newVertex
            case 13 =>
              val newVertex = Node(x, y, topVertices(x))
              leftVertex = newVertex

            case 14 =>
              val newVertex = Node(x, y, leftVertex)
              topVertices(x).next = newVertex


          }
          offset += 1
        }
        offset = offset - (resolution + 2 * bufferSize + 1) + resolution
      }
      //ring formation
      val factory: GeometryFactory = GeometryReader.DefaultGeometryFactory
      val rings = new ArrayBuffer[LinearRing]()
      for (corner <- corners; if !corner.visited) {
        //  First, count the number of corners to prepare a CoordinateSequence of the right size
        var iCorner = 0
        var p = corner
        do {
          p = p.next
          iCorner += 1
        } while (p != corner)
        val coords = factory.getCoordinateSequenceFactory.create(iCorner + 1, 2)
        p = corner
        iCorner = 0
        do {
          coords.setOrdinate(iCorner, 0, p.x)
          coords.setOrdinate(iCorner, 1, p.y)
          p.visited = true
          p = p.next
          iCorner += 1
        } while (p != corner)
        // Make last coordinate similar to the first one
        coords.setOrdinate(iCorner, 0, p.x)
        coords.setOrdinate(iCorner, 1, p.y)
        // Create the linear ring
        rings.append(factory.createLinearRing(coords))

      }
      rings.toArray
      val createdPolygons = ArrayBuffer[Polygon]()
      for(ring<- rings){
        createdPolygons.append(factory.createPolygon(ring,null))
      }
      val polygons=factory.createMultiPolygon(createdPolygons.toArray)
      vectorLayerBuilder.addFeature(Feature.create(null, polygons))
      vectorLayerBuilder.build()
    }
    val vectorTileBuilder = VectorTile.Tile.newBuilder()
    vectorTileBuilder.addLayers(layer)
    vectorTileBuilder.build()
  }
}



object IntermediateTile {
  /** The maximum number of points in one vector tile before we decide to rasterize */
  //val maxPointsPerTile: Long = 1000000000000L


  val maxPointsPerTile: Long = 50L
}
case class Node(x: Int, y: Int, var next: Node, var visited: Boolean = false)